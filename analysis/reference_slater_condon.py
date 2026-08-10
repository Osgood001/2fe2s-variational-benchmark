"""Small, explicit Slater--Condon reference implementation.

This module deliberately does not depend on PySCF.  It parses the real
FCIDUMP, validates a sparse fixed-particle wave function, and evaluates
Hamiltonian matrix elements by applying fermionic creation and annihilation
operators to spin-orbital bit strings.  It is used as an independent numerical
path for the 512-determinant variational control.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


NORB = 20
NALPHA = 15
NBETA = 15


class InvalidWavefunction(ValueError):
    """Raised when the sparse state does not satisfy the declared sector."""


def _eri_key(i: int, j: int, k: int, l: int):
    pair_a = tuple(sorted((i, j), reverse=True))
    pair_b = tuple(sorted((k, l), reverse=True))
    return tuple(sorted((pair_a, pair_b), reverse=True))


def parse_fcidump(path: Path):
    """Parse one- and two-electron integrals in eightfold FCIDUMP symmetry."""
    one = {}
    two = {}
    core = 0.0
    in_header = True
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if in_header:
            if line.upper().startswith("&END") or line == "/":
                in_header = False
            continue
        fields = line.replace("D", "E").split()
        if len(fields) != 5:
            continue
        value = float(fields[0])
        i, j, k, l = (int(token) for token in fields[1:])
        if i == j == k == l == 0:
            core = value
        elif k == l == 0:
            one[tuple(sorted((i - 1, j - 1), reverse=True))] = value
        else:
            two[_eri_key(i - 1, j - 1, k - 1, l - 1)] = value
    return one, two, core


def _h(one, p: int, q: int) -> float:
    return one.get(tuple(sorted((p, q), reverse=True)), 0.0)


def _eri(two, p: int, q: int, r: int, s: int) -> float:
    return two.get(_eri_key(p, q, r, s), 0.0)


def _annihilate(bits: int, orbital: int):
    mask = 1 << orbital
    if not bits & mask:
        return None
    phase = -1.0 if (bits & (mask - 1)).bit_count() % 2 else 1.0
    return bits ^ mask, phase


def _create(bits: int, orbital: int):
    mask = 1 << orbital
    if bits & mask:
        return None
    phase = -1.0 if (bits & (mask - 1)).bit_count() % 2 else 1.0
    return bits | mask, phase


def _spin_orbitals(bits: int):
    return [index for index in range(2 * NORB) if bits >> index & 1]


def _one_term(one, bra: int, ket: int, p: int, q: int) -> float:
    if p // NORB != q // NORB:
        return 0.0
    result = _annihilate(ket, q)
    if result is None:
        return 0.0
    state, phase = result
    result = _create(state, p)
    if result is None or result[0] != bra:
        return 0.0
    return phase * result[1] * _h(one, p % NORB, q % NORB)


def _two_term(two, bra: int, ket: int, p: int, q: int, r: int, s: int) -> float:
    if p // NORB != q // NORB or r // NORB != s // NORB:
        return 0.0
    state = ket
    phase = 1.0
    for operation, orbital in (("a", q), ("a", s), ("c", r), ("c", p)):
        result = (
            _annihilate(state, orbital)
            if operation == "a"
            else _create(state, orbital)
        )
        if result is None:
            return 0.0
        state, local_phase = result
        phase *= local_phase
    if state != bra:
        return 0.0
    return 0.5 * phase * _eri(two, p % NORB, q % NORB, r % NORB, s % NORB)


def matrix_element(one, two, core: float, bra: int, ket: int) -> float:
    """Return <bra|H|ket> by explicit second-quantized operator action."""
    removed = _spin_orbitals(ket & ~bra)
    added = _spin_orbitals(bra & ~ket)
    if len(removed) != len(added) or len(removed) > 2:
        return 0.0
    common = _spin_orbitals(bra & ket)
    value = core if not removed else 0.0
    if not removed:
        for q in common:
            value += _one_term(one, bra, ket, q, q)
        for q in common:
            for s in common:
                if q == s:
                    continue
                value += _two_term(two, bra, ket, q, q, s, s)
                value += _two_term(two, bra, ket, s, q, q, s)
    elif len(removed) == 1:
        q = removed[0]
        p = added[0]
        value += _one_term(one, bra, ket, p, q)
        for orbital in common:
            for annihilated_q, annihilated_s in ((q, orbital), (orbital, q)):
                for created_p, created_r in ((p, orbital), (orbital, p)):
                    value += _two_term(
                        two,
                        bra,
                        ket,
                        created_p,
                        annihilated_q,
                        created_r,
                        annihilated_s,
                    )
    else:
        for q, s in ((removed[0], removed[1]), (removed[1], removed[0])):
            for p, r in ((added[0], added[1]), (added[1], added[0])):
                value += _two_term(two, bra, ket, p, q, r, s)
    return value


def _parse_bits(text: str) -> int:
    value = text.strip()
    if len(value) != NORB or any(character not in "01" for character in value):
        raise InvalidWavefunction("occupations must be 20-character bit strings")
    return sum((character == "1") << index for index, character in enumerate(value))


def load_wavefunction(path: Path, expected_determinants: int = 512):
    """Load and validate alpha/beta occupation strings and real coefficients."""
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["alpha", "beta", "coefficient"]:
            raise InvalidWavefunction("expected alpha,beta,coefficient columns")
        for row in reader:
            alpha = _parse_bits(row["alpha"])
            beta = _parse_bits(row["beta"])
            if alpha.bit_count() != NALPHA or beta.bit_count() != NBETA:
                raise InvalidWavefunction("determinant lies outside the (15,15) sector")
            coefficient = float(row["coefficient"])
            if not math.isfinite(coefficient):
                raise InvalidWavefunction("coefficient is not finite")
            rows.append((alpha | (beta << NORB), coefficient))
    if len(rows) != expected_determinants:
        raise InvalidWavefunction(
            f"expected {expected_determinants} determinants, found {len(rows)}"
        )
    if len({bits for bits, _ in rows}) != len(rows):
        raise InvalidWavefunction("duplicate determinant")
    norm = sum(coefficient * coefficient for _, coefficient in rows)
    if not math.isfinite(norm) or norm < 1e-16:
        raise InvalidWavefunction("wave-function norm is zero or invalid")
    return rows, norm
