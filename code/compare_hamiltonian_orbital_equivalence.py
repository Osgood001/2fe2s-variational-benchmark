#!/usr/bin/env python3
"""Prove orbital-equivalence of two real FCIDUMP Hamiltonians.

The one-electron matrices have nondegenerate spectra for this instance.  Their
ordered eigenvectors therefore define a common canonical orbital basis up to
independent column signs.  We infer those signs from two-electron integrals and
then compare the entire transformed rank-four tensor.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from pyscf import ao2mo
from pyscf.tools import fcidump


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: compare_hamiltonian_orbital_equivalence.py A B OUT")
    path_a, path_b, output_path = map(Path, sys.argv[1:])
    data_a = fcidump.read(str(path_a), verbose=False)
    data_b = fcidump.read(str(path_b), verbose=False)

    headers = {}
    for key in ("NORB", "NELEC", "MS2", "ECORE"):
        value_a = data_a[key]
        value_b = data_b[key]
        headers[key] = {
            "a": value_a,
            "b": value_b,
            "equal": bool(value_a == value_b),
        }
    if not all(record["equal"] for record in headers.values()):
        raise RuntimeError("FCIDUMP sectors or constant energies differ")

    norb = int(data_a["NORB"])
    h1_a = np.asarray(data_a["H1"], dtype=np.float64)
    h1_b = np.asarray(data_b["H1"], dtype=np.float64)
    eri_a = ao2mo.restore(1, np.asarray(data_a["H2"]), norb)
    eri_b = ao2mo.restore(1, np.asarray(data_b["H2"]), norb)

    eigenvalues_a, eigenvectors_a = np.linalg.eigh(h1_a)
    eigenvalues_b, eigenvectors_b = np.linalg.eigh(h1_b)
    minimum_h1_gap = float(np.min(np.diff(eigenvalues_a)))
    if minimum_h1_gap < 1e-8:
        raise RuntimeError("one-electron spectrum is too degenerate for this proof")

    canonical_a = np.einsum(
        "pi,qj,rk,sl,pqrs->ijkl",
        eigenvectors_a,
        eigenvectors_a,
        eigenvectors_a,
        eigenvectors_a,
        eri_a,
        optimize=True,
    )
    canonical_b = np.einsum(
        "pi,qj,rk,sl,pqrs->ijkl",
        eigenvectors_b,
        eigenvectors_b,
        eigenvectors_b,
        eigenvectors_b,
        eri_b,
        optimize=True,
    )

    # Eigenvector columns are individually sign-indeterminate.  Fix a gauge
    # s[0]=+1, then determine each remaining relative sign from the largest
    # available integral g[i,r,r,r] connected to an already resolved r.
    signs = np.ones(norb, dtype=np.float64)
    known = {0}
    sign_witnesses = []
    while len(known) < norb:
        best = None
        for i in set(range(norb)) - known:
            for reference in known:
                magnitude = min(
                    abs(canonical_a[i, reference, reference, reference]),
                    abs(canonical_b[i, reference, reference, reference]),
                )
                if best is None or magnitude > best[0]:
                    best = (magnitude, i, reference)
        magnitude, i, reference = best
        if magnitude < 1e-12:
            raise RuntimeError("could not determine all eigenvector signs")
        ratio = (
            canonical_a[i, reference, reference, reference]
            / canonical_b[i, reference, reference, reference]
        )
        signs[i] = np.sign(ratio) * signs[reference]
        known.add(i)
        sign_witnesses.append(
            {
                "orbital": i,
                "reference": reference,
                "minimum_abs_integral": float(magnitude),
            }
        )

    sign_product = np.einsum("i,j,k,l->ijkl", signs, signs, signs, signs)
    tensor_difference = canonical_a - sign_product * canonical_b
    orbital_rotation_b_to_a = eigenvectors_a @ np.diag(signs) @ eigenvectors_b.T
    rotated_h1_b = orbital_rotation_b_to_a @ h1_b @ orbital_rotation_b_to_a.T

    result = {
        "schema_version": "fcidump-orbital-equivalence/v1",
        "path_a": str(path_a.resolve()),
        "path_b": str(path_b.resolve()),
        "sha256_a": sha256_file(path_a),
        "sha256_b": sha256_file(path_b),
        "byte_identical": bool(path_a.read_bytes() == path_b.read_bytes()),
        "headers": headers,
        "minimum_h1_eigenvalue_gap": minimum_h1_gap,
        "h1_eigenvalue_max_abs_difference": float(
            np.max(np.abs(eigenvalues_a - eigenvalues_b))
        ),
        "h1_rotation_max_abs_difference": float(
            np.max(np.abs(h1_a - rotated_h1_b))
        ),
        "eri_frobenius_a": float(np.linalg.norm(eri_a)),
        "eri_frobenius_b": float(np.linalg.norm(eri_b)),
        "canonical_eri_max_abs_difference": float(
            np.max(np.abs(tensor_difference))
        ),
        "canonical_eri_rms_difference": float(
            np.sqrt(np.mean(tensor_difference**2))
        ),
        "canonical_eri_allclose_rtol0_atol1e-10": bool(
            np.allclose(canonical_a, sign_product * canonical_b, rtol=0, atol=1e-10)
        ),
        "rotation_orthogonality_max_abs_difference": float(
            np.max(np.abs(orbital_rotation_b_to_a.T @ orbital_rotation_b_to_a - np.eye(norb)))
        ),
        "eigenvector_signs": signs.astype(int).tolist(),
        "sign_witnesses": sign_witnesses,
        "interpretation": (
            "The files differ in orbital basis but represent the same many-electron "
            "Hamiltonian to FCIDUMP numerical precision."
        ),
    }
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
