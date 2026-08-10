#!/usr/bin/env python3
"""Transform a fixed-electron CI tensor to a rotated orbital basis exactly.

For 15 electrons in 20 orbitals, direct determinant transformation requires
15x15 minors.  Jacobi's complementary-minor identity maps the problem to the
five-hole representation, where PySCF evaluates only 5x5 minors.  The same
physical Ritz state is then contracted with the independently transformed
FCIDUMP, providing an energy/residual invariance check before any top-K study.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import pyscf
from pyscf import ao2mo, fci, lib
from pyscf.tools import fcidump


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def load_rotation(path: Path, norb: int) -> np.ndarray:
    rotation = np.empty((norb, norb), dtype=np.float64)
    seen: set[tuple[int, int]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["source_orbital", "natural_orbital", "coefficient"]
        if reader.fieldnames != expected:
            raise SystemExit(f"rotation columns mismatch: {reader.fieldnames}")
        for row in reader:
            source = int(row["source_orbital"]) - 1
            target = int(row["natural_orbital"]) - 1
            if not 0 <= source < norb or not 0 <= target < norb:
                raise SystemExit("rotation index out of range")
            if (source, target) in seen:
                raise SystemExit("duplicate rotation entry")
            seen.add((source, target))
            rotation[source, target] = float(row["coefficient"])
    if len(seen) != norb * norb:
        raise SystemExit("rotation matrix is incomplete")
    return rotation


def atomic_write_json(path: Path, value: dict[str, object]) -> None:
    """Write a diagnostic atomically so failed scientific gates remain inspectable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def address_maps(norb: int, nelec: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return electron->hole address, its inverse, and electron parity signs."""
    nhole = norb - nelec
    electron_strings = np.asarray(
        fci.cistring.make_strings(range(norb), nelec), dtype=np.int64
    )
    full_mask = (1 << norb) - 1
    hole_strings = np.bitwise_xor(electron_strings, full_mask)
    hole_for_electron = np.asarray(
        fci.cistring.strs2addr(norb, nhole, hole_strings), dtype=np.int64
    )
    electron_for_hole = np.empty_like(hole_for_electron)
    electron_for_hole[hole_for_electron] = np.arange(electron_strings.size)
    signs = np.empty(electron_strings.size, dtype=np.float64)
    for address, string in enumerate(electron_strings):
        orbital_sum = sum(index for index in range(norb) if int(string) >> index & 1)
        signs[address] = -1.0 if orbital_sum % 2 else 1.0
    return hole_for_electron, electron_for_hole, signs


def transform_via_holes(
    electron_ci: np.ndarray,
    norb: int,
    nelec: tuple[int, int],
    rotation: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    if nelec[0] != nelec[1]:
        raise SystemExit("this optimized path requires equal alpha/beta populations")
    if electron_ci.shape[0] != electron_ci.shape[1]:
        raise SystemExit("this optimized path requires equal spin-string dimensions")
    hole_for_electron, electron_for_hole, electron_sign = address_maps(
        norb, nelec[0]
    )
    sign_for_hole = electron_sign[electron_for_hole]

    mapping_started = time.time()
    hole_ci = electron_ci[np.ix_(electron_for_hole, electron_for_hole)]
    hole_ci *= sign_for_hole[:, None]
    hole_ci *= sign_for_hole[None, :]
    mapping_to_holes_s = time.time() - mapping_started

    transform_started = time.time()
    transformed_hole_ci = fci.addons.transform_ci(
        hole_ci, (norb - nelec[0], norb - nelec[1]), rotation
    )
    transform_s = time.time() - transform_started
    del hole_ci
    gc.collect()

    mapping_started = time.time()
    transformed_electron_ci = transformed_hole_ci[
        np.ix_(hole_for_electron, hole_for_electron)
    ]
    transformed_electron_ci *= electron_sign[:, None]
    transformed_electron_ci *= electron_sign[None, :]
    transformed_electron_ci *= float(np.linalg.det(rotation) ** 2)
    mapping_from_holes_s = time.time() - mapping_started
    del transformed_hole_ci
    gc.collect()
    return transformed_electron_ci, {
        "mapping_to_holes_s": mapping_to_holes_s,
        "five_hole_transform_s": transform_s,
        "mapping_from_holes_s": mapping_from_holes_s,
    }


def small_system_gate() -> dict[str, float | bool]:
    rng = np.random.default_rng(20260809)
    norb = 6
    nelec = (4, 4)
    random_matrix = rng.normal(size=(norb, norb))
    rotation, _ = np.linalg.qr(random_matrix)
    shape = (
        fci.cistring.num_strings(norb, nelec[0]),
        fci.cistring.num_strings(norb, nelec[1]),
    )
    ci = rng.normal(size=shape)
    ci /= np.linalg.norm(ci)
    direct = fci.addons.transform_ci(ci, nelec, rotation)
    complementary, _ = transform_via_holes(ci, norb, nelec, rotation)
    max_abs = float(np.max(np.abs(direct - complementary)))
    return {
        "norb": norb,
        "nelec_alpha": nelec[0],
        "max_abs_direct_vs_complementary_minor": max_abs,
        "passed": max_abs <= 2e-12,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_fcidump", type=Path)
    parser.add_argument("natural_fcidump", type=Path)
    parser.add_argument("source_vector", type=Path)
    parser.add_argument("rotation_csv", type=Path)
    parser.add_argument("output_vector", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--expected-source-fcidump-sha256", required=True)
    parser.add_argument("--expected-natural-fcidump-sha256", required=True)
    parser.add_argument("--expected-source-vector-sha256", required=True)
    parser.add_argument("--expected-rotation-sha256", required=True)
    parser.add_argument("--reference-energy-eh", type=float, required=True)
    parser.add_argument("--reference-residual-eh", type=float, required=True)
    args = parser.parse_args()
    started = time.time()
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "64")))

    hashes = {
        "source_fcidump_sha256": sha256_file(args.source_fcidump),
        "natural_fcidump_sha256": sha256_file(args.natural_fcidump),
        "source_vector_sha256": sha256_file(args.source_vector),
        "rotation_sha256": sha256_file(args.rotation_csv),
    }
    expected_hashes = {
        "source_fcidump_sha256": args.expected_source_fcidump_sha256,
        "natural_fcidump_sha256": args.expected_natural_fcidump_sha256,
        "source_vector_sha256": args.expected_source_vector_sha256,
        "rotation_sha256": args.expected_rotation_sha256,
    }
    if hashes != expected_hashes:
        raise SystemExit(f"input hash gate failed: {hashes} != {expected_hashes}")

    source_data = fcidump.read(str(args.source_fcidump), verbose=False)
    natural_data = fcidump.read(str(args.natural_fcidump), verbose=False)
    norb = int(source_data["NORB"])
    nelec_total = int(source_data["NELEC"])
    ms2 = int(source_data["MS2"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    if (
        int(natural_data["NORB"]) != norb
        or int(natural_data["NELEC"]) != nelec_total
        or int(natural_data["MS2"]) != ms2
    ):
        raise SystemExit("source and natural-basis FCIDUMP sectors differ")
    rotation = load_rotation(args.rotation_csv, norb)
    orthogonality_error = float(
        np.max(np.abs(rotation.T @ rotation - np.eye(norb)))
    )
    if orthogonality_error > 2e-12:
        raise SystemExit(f"rotation orthogonality gate failed: {orthogonality_error}")
    small_gate = small_system_gate()
    if not small_gate["passed"]:
        raise SystemExit("complementary-minor small-system gate failed")

    source_vector = np.load(args.source_vector, allow_pickle=False)
    expected_shape = (
        fci.cistring.num_strings(norb, nelec[0]),
        fci.cistring.num_strings(norb, nelec[1]),
    )
    if source_vector.shape != expected_shape or source_vector.dtype != np.float64:
        raise SystemExit("source vector shape/dtype mismatch")
    source_norm2 = float(np.dot(source_vector.reshape(-1), source_vector.reshape(-1)))
    transformed_vector, timing = transform_via_holes(
        source_vector, norb, nelec, rotation
    )
    transformed_norm2 = float(
        np.dot(transformed_vector.reshape(-1), transformed_vector.reshape(-1))
    )
    args.output_vector.parent.mkdir(parents=True, exist_ok=True)
    np.save(args.output_vector, transformed_vector, allow_pickle=False)
    transformed_vector_sha256 = sha256_file(args.output_vector)

    # Contract the source representation in this same hash-locked process.
    # This makes the invariance test self-contained rather than comparing the
    # transformed result only with command-line reference constants.
    source_h1e = np.asarray(source_data["H1"], dtype=np.float64)
    source_eri = ao2mo.restore(
        8, np.asarray(source_data["H2"], dtype=np.float64), norb
    )
    source_ecore = float(source_data["ECORE"])
    source_h2e = fci.direct_spin1.absorb_h1e(
        source_h1e, source_eri, norb, nelec, 0.5
    )
    source_contracted = np.asarray(
        fci.direct_spin1.contract_2e(
            source_h2e, source_vector, norb, nelec
        ),
        dtype=np.float64,
    ).reshape(expected_shape)
    source_electronic_energy = float(
        np.dot(source_vector.reshape(-1), source_contracted.reshape(-1))
        / source_norm2
    )
    source_energy = source_electronic_energy + source_ecore
    source_residual = float(
        np.linalg.norm(
            (
                source_contracted
                - source_electronic_energy * source_vector
            ).reshape(-1)
        )
        / np.sqrt(source_norm2)
    )
    del source_contracted, source_h2e, source_h1e, source_eri, source_data
    gc.collect()

    natural_h1e = np.asarray(natural_data["H1"], dtype=np.float64)
    natural_eri = ao2mo.restore(
        8, np.asarray(natural_data["H2"], dtype=np.float64), norb
    )
    natural_ecore = float(natural_data["ECORE"])
    natural_h2e = fci.direct_spin1.absorb_h1e(
        natural_h1e, natural_eri, norb, nelec, 0.5
    )
    contracted = np.asarray(
        fci.direct_spin1.contract_2e(
            natural_h2e, transformed_vector, norb, nelec
        ),
        dtype=np.float64,
    ).reshape(expected_shape)
    electronic_energy = float(
        np.dot(transformed_vector.reshape(-1), contracted.reshape(-1))
        / transformed_norm2
    )
    energy = electronic_energy + natural_ecore
    residual = float(
        np.linalg.norm(
            (contracted - electronic_energy * transformed_vector).reshape(-1)
        )
        / np.sqrt(transformed_norm2)
    )
    gates = {
        "small_system_complementary_minor": bool(small_gate["passed"]),
        "rotation_orthogonality": orthogonality_error <= 2e-12,
        "norm_invariance": abs(transformed_norm2 - source_norm2) <= 2e-10,
        "source_reference_energy": abs(
            source_energy - args.reference_energy_eh
        ) <= 2e-9,
        "source_reference_residual": abs(
            source_residual - args.reference_residual_eh
        ) <= 2e-8,
        "same_job_energy_invariance": abs(energy - source_energy) <= 2e-9,
        "same_job_residual_invariance": abs(residual - source_residual) <= 2e-8,
        "transformed_reference_energy": abs(
            energy - args.reference_energy_eh
        ) <= 2e-9,
        "transformed_reference_residual": abs(
            residual - args.reference_residual_eh
        ) <= 2e-8,
    }
    result = {
        "schema_version": "2fe2s-ci-orbital-transform/v2",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": pyscf.__version__,
        "threads": lib.num_threads(),
        **hashes,
        "output_vector_sha256": transformed_vector_sha256,
        "output_vector_bytes": args.output_vector.stat().st_size,
        "output_vector_shape": list(transformed_vector.shape),
        "rotation_orthogonality_max_abs": orthogonality_error,
        "rotation_determinant": float(np.linalg.det(rotation)),
        "source_norm2": source_norm2,
        "transformed_norm2": transformed_norm2,
        "source_energy_eh": source_energy,
        "source_energy_difference_from_reference_eh": source_energy
        - args.reference_energy_eh,
        "source_residual_eh": source_residual,
        "source_residual_difference_from_reference_eh": source_residual
        - args.reference_residual_eh,
        "transformed_energy_eh": energy,
        "energy_difference_natural_minus_source_eh": energy - source_energy,
        "energy_difference_from_reference_eh": energy
        - args.reference_energy_eh,
        "transformed_residual_eh": residual,
        "residual_difference_natural_minus_source_eh": residual
        - source_residual,
        "residual_difference_from_reference_eh": residual
        - args.reference_residual_eh,
        "small_system_gate": small_gate,
        "timing": timing,
        "wall_s": time.time() - started,
        "gates": gates,
    }
    atomic_write_json(args.output_json, result)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    if not all(gates.values()):
        raise SystemExit("orbital-transform gate failed: " + json.dumps(gates))


if __name__ == "__main__":
    main()
