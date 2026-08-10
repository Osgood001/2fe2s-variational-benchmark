#!/usr/bin/env python3
"""Independent validation of a Hamiltonian-optimized 512-determinant state."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fcidump", type=Path)
    parser.add_argument("wavefunction", type=Path)
    parser.add_argument("reference_implementation", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.time()
    lib.num_threads(64)

    spec = importlib.util.spec_from_file_location(
        "reference_slater_condon", args.reference_implementation
    )
    reference = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(reference)
    rows, clean_norm = reference.load_wavefunction(args.wavefunction)
    one, two, core = reference.parse_fcidump(args.fcidump)
    clean_numerator = 0.0
    for index, (bra, coefficient_i) in enumerate(rows):
        clean_numerator += (
            coefficient_i
            * coefficient_i
            * reference.matrix_element(one, two, core, bra, bra)
        )
        for ket, coefficient_j in rows[:index]:
            clean_numerator += (
                2.0
                * coefficient_i
                * coefficient_j
                * reference.matrix_element(one, two, core, bra, ket)
            )
    clean_energy = clean_numerator / clean_norm

    data = fcidump.read(str(args.fcidump), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    ecore = float(data["ECORE"])
    h1e = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)
    nstr_alpha = fci.cistring.num_strings(norb, nelec[0])
    nstr_beta = fci.cistring.num_strings(norb, nelec[1])
    vector = np.zeros((nstr_alpha, nstr_beta), dtype=np.float64)
    support: list[tuple[int, int]] = []
    coefficients: list[float] = []
    orbital_mask = (1 << norb) - 1
    for determinant, coefficient in rows:
        alpha = determinant & orbital_mask
        beta = determinant >> norb
        alpha_address = int(fci.cistring.str2addr(norb, nelec[0], alpha))
        beta_address = int(fci.cistring.str2addr(norb, nelec[1], beta))
        vector[alpha_address, beta_address] = coefficient
        support.append((alpha_address, beta_address))
        coefficients.append(coefficient)

    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, nelec, 0.5)
    contracted = np.asarray(
        fci.direct_spin1.contract_2e(h2e, vector, norb, nelec), dtype=np.float64
    ).reshape(vector.shape)
    norm = float(np.dot(vector.reshape(-1), vector.reshape(-1)))
    electronic_energy = float(
        np.dot(vector.reshape(-1), contracted.reshape(-1)) / norm
    )
    pyscf_energy = electronic_energy + ecore
    full_residual = float(
        np.linalg.norm((contracted - electronic_energy * vector).reshape(-1))
        / np.sqrt(norm)
    )
    support_hc = np.asarray(
        [contracted[alpha_address, beta_address] for alpha_address, beta_address in support]
    )
    coefficient_array = np.asarray(coefficients, dtype=np.float64)
    projected_residual = float(
        np.linalg.norm(support_hc - electronic_energy * coefficient_array)
        / np.linalg.norm(coefficient_array)
    )
    result = {
        "schema_version": "2fe2s-optimized-512-validation/v1",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": __import__("pyscf").__version__,
        "threads": lib.num_threads(),
        "fcidump_sha256": sha256_file(args.fcidump),
        "wavefunction_sha256": sha256_file(args.wavefunction),
        "wavefunction_bytes": args.wavefunction.stat().st_size,
        "determinants": len(rows),
        "norm2_clean_room": clean_norm,
        "norm2_pyscf_embedding": norm,
        "energy_clean_room_slater_condon_eh": clean_energy,
        "energy_pyscf_direct_contraction_eh": pyscf_energy,
        "energy_path_difference_eh": pyscf_energy - clean_energy,
        "projected_subspace_residual_eh": projected_residual,
        "full_space_residual_eh": full_residual,
        "wall_s": time.time() - started,
        "interpretation": (
            "The projected residual tests coefficient stationarity in the fixed "
            "512-determinant support; the full residual is not expected to vanish."
        ),
        "gates": {
            "determinant_count": len(rows) == 512,
            "norm": abs(norm - 1.0) <= 2e-12 and abs(clean_norm - 1.0) <= 2e-12,
            "independent_energy_agreement": abs(pyscf_energy - clean_energy) <= 1e-10,
            "projected_stationarity": projected_residual <= 1e-9,
        },
    }
    if not all(result["gates"].values()):
        raise SystemExit("validation gate failed: " + json.dumps(result["gates"]))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
