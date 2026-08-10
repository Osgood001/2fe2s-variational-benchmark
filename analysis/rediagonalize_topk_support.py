#!/usr/bin/env python3
"""Rediagonalize the amplitude-ranked top-K support of a full CI state.

The support is selected only from the hash-locked full-state amplitudes.  Its
Hamiltonian matrix is then built by the explicit fermionic-operator reference
implementation and diagonalized exactly.  A separate PySCF full-space
contraction validates the optimized coefficients and projected stationarity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
from pyscf import ao2mo, fci, lib
from pyscf.tools import fcidump


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def occupation_string(bits: int, norb: int) -> str:
    return "".join("1" if bits >> orbital & 1 else "0" for orbital in range(norb))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fcidump", type=Path)
    parser.add_argument("full_vector", type=Path)
    parser.add_argument("reference_implementation", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--k", type=int, default=512)
    parser.add_argument("--expected-fcidump-sha256", required=True)
    parser.add_argument("--expected-vector-sha256", required=True)
    parser.add_argument("--expected-reference-sha256", required=True)
    parser.add_argument("--reference-topk-weight", type=float, required=True)
    parser.add_argument("--reference-fixed-energy-eh", type=float, required=True)
    args = parser.parse_args()
    started = time.time()
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "64")))

    hashes = {
        "fcidump_sha256": sha256_file(args.fcidump),
        "full_vector_sha256": sha256_file(args.full_vector),
        "reference_implementation_sha256": sha256_file(
            args.reference_implementation
        ),
    }
    expected = {
        "fcidump_sha256": args.expected_fcidump_sha256,
        "full_vector_sha256": args.expected_vector_sha256,
        "reference_implementation_sha256": args.expected_reference_sha256,
    }
    if hashes != expected:
        raise SystemExit(f"input hash gate failed: {hashes} != {expected}")

    spec = importlib.util.spec_from_file_location(
        "reference_slater_condon", args.reference_implementation
    )
    reference = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(reference)
    one, two, core = reference.parse_fcidump(args.fcidump)

    data = fcidump.read(str(args.fcidump), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    ecore = float(data["ECORE"])
    nalpha = fci.cistring.num_strings(norb, nelec[0])
    nbeta = fci.cistring.num_strings(norb, nelec[1])
    vector = np.load(args.full_vector, allow_pickle=False)
    if vector.shape != (nalpha, nbeta) or vector.dtype != np.float64:
        raise SystemExit("full-vector shape/dtype gate failed")
    flat = vector.reshape(-1)

    weights = np.square(flat, dtype=np.float64)
    selected = np.argpartition(weights, weights.size - args.k)[-args.k:]
    selected = np.sort(selected.astype(np.int64, copy=False))
    retained_weight = float(np.sum(weights[selected], dtype=np.float64))
    fixed_coefficients = np.asarray(flat[selected], dtype=np.float64)
    fixed_coefficients /= np.linalg.norm(fixed_coefficients)
    del weights

    alpha_strings = np.asarray(
        fci.cistring.make_strings(range(norb), nelec[0]), dtype=np.int64
    )
    beta_strings = np.asarray(
        fci.cistring.make_strings(range(norb), nelec[1]), dtype=np.int64
    )
    alpha_addresses = selected // nbeta
    beta_addresses = selected % nbeta
    determinant_bits = [
        int(alpha_strings[alpha_address])
        | (int(beta_strings[beta_address]) << norb)
        for alpha_address, beta_address in zip(alpha_addresses, beta_addresses)
    ]

    hamiltonian = np.empty((args.k, args.k), dtype=np.float64)
    for row, bra in enumerate(determinant_bits):
        hamiltonian[row, row] = reference.matrix_element(
            one, two, core, bra, bra
        )
        for column in range(row):
            value = reference.matrix_element(
                one, two, core, bra, determinant_bits[column]
            )
            hamiltonian[row, column] = value
            hamiltonian[column, row] = value
        if (row + 1) % 64 == 0 or row + 1 == args.k:
            print(f"REFERENCE_MATRIX rows={row + 1}/{args.k}", flush=True)

    fixed_energy = float(
        fixed_coefficients @ hamiltonian @ fixed_coefficients
    )
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    optimized_energy_reference = float(eigenvalues[0])
    optimized_coefficients = np.asarray(eigenvectors[:, 0], dtype=np.float64)
    pivot = int(np.argmax(np.abs(optimized_coefficients)))
    if optimized_coefficients[pivot] < 0.0:
        optimized_coefficients *= -1.0
    projected_residual_reference = float(
        np.linalg.norm(
            hamiltonian @ optimized_coefficients
            - optimized_energy_reference * optimized_coefficients
        )
    )
    fixed_optimized_overlap = float(
        abs(np.dot(fixed_coefficients, optimized_coefficients))
    )

    sparse_vector = np.zeros_like(vector)
    sparse_vector.reshape(-1)[selected] = optimized_coefficients
    h1e = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)
    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, nelec, 0.5)
    contracted = np.asarray(
        fci.direct_spin1.contract_2e(h2e, sparse_vector, norb, nelec),
        dtype=np.float64,
    ).reshape(sparse_vector.shape)
    norm2 = float(np.dot(optimized_coefficients, optimized_coefficients))
    electronic_energy = float(
        np.dot(sparse_vector.reshape(-1), contracted.reshape(-1)) / norm2
    )
    optimized_energy_pyscf = electronic_energy + ecore
    full_space_residual = float(
        np.linalg.norm(
            (contracted - electronic_energy * sparse_vector).reshape(-1)
        )
        / np.sqrt(norm2)
    )
    support_hc_total = contracted.reshape(-1)[selected] + ecore * optimized_coefficients
    projected_residual_pyscf = float(
        np.linalg.norm(
            support_hc_total
            - optimized_energy_pyscf * optimized_coefficients
        )
        / np.sqrt(norm2)
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["alpha", "beta", "coefficient"])
        for alpha_address, beta_address, coefficient in zip(
            alpha_addresses, beta_addresses, optimized_coefficients
        ):
            writer.writerow(
                [
                    occupation_string(int(alpha_strings[alpha_address]), norb),
                    occupation_string(int(beta_strings[beta_address]), norb),
                    format(float(coefficient), ".17g"),
                ]
            )
    output_csv_sha256 = sha256_file(args.output_csv)

    gates = {
        "topk_weight_reproduced": abs(
            retained_weight - args.reference_topk_weight
        ) <= 2e-12,
        "fixed_energy_reproduced": abs(
            fixed_energy - args.reference_fixed_energy_eh
        ) <= 1e-10,
        "reference_projected_stationarity": projected_residual_reference
        <= 1e-10,
        "pyscf_projected_stationarity": projected_residual_pyscf <= 1e-9,
        "independent_energy_agreement": abs(
            optimized_energy_pyscf - optimized_energy_reference
        ) <= 1e-10,
        "coefficient_norm": abs(norm2 - 1.0) <= 2e-12,
        "rediagonalization_is_variational": optimized_energy_reference
        <= fixed_energy + 1e-12,
    }
    result = {
        "schema_version": "2fe2s-topk-rediagonalization/v1",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": __import__("pyscf").__version__,
        "threads": lib.num_threads(),
        **hashes,
        "k": args.k,
        "retained_norm2": retained_weight,
        "fixed_coefficient_energy_reference_eh": fixed_energy,
        "same_support_optimized_energy_reference_eh": optimized_energy_reference,
        "same_support_optimized_energy_pyscf_eh": optimized_energy_pyscf,
        "energy_gain_from_coefficient_reoptimization_eh": fixed_energy
        - optimized_energy_reference,
        "fixed_vs_optimized_absolute_overlap": fixed_optimized_overlap,
        "projected_residual_reference_eh": projected_residual_reference,
        "projected_residual_pyscf_eh": projected_residual_pyscf,
        "full_space_residual_pyscf_eh": full_space_residual,
        "optimized_norm2": norm2,
        "output_wavefunction_sha256": output_csv_sha256,
        "output_wavefunction_bytes": args.output_csv.stat().st_size,
        "wall_s": time.time() - started,
        "interpretation": (
            "Only coefficients were optimized; the determinant support is the "
            "amplitude-ranked top-K support of the frozen full Ritz vector."
        ),
        "gates": gates,
    }
    atomic_write_json(args.output_json, result)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    if not all(gates.values()):
        raise SystemExit("rediagonalization gate failed: " + json.dumps(gates))


if __name__ == "__main__":
    main()
