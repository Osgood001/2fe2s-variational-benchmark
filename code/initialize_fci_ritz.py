#!/usr/bin/env python3
"""Reproduce and preserve a variational FCI Davidson Ritz vector for 2Fe-2S.

This intentionally stops after a declared number of Davidson cycles.  A stopped
Ritz vector can still define a variational upper bound, but it must not be called
"converged" unless PySCF says so.  The full vector, its hash, its directly
contracted Rayleigh quotient, and its residual are preserved for reproduction.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import pyscf
from pyscf import ao2mo, fci, lib
from pyscf.tools import fcidump


PUBLIC_HCI_BOUND_EH = -116.605425
EXPECTED_FCIDUMP_SHA256 = os.environ.get(
    "EXPECTED_FCIDUMP_SHA256",
    "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7",
)


def sha256_file(path: Path, block_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        block = handle.read(block_size)
        while block:
            digest.update(block)
            block = handle.read(block_size)
    return digest.hexdigest()


def emit(event: str, **values: object) -> None:
    payload = {"elapsed_s": round(time.time() - START, 3), "event": event, **values}
    print("SIMULATION_JSON " + json.dumps(payload, sort_keys=True), flush=True)


START = time.time()


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: initialize_fci_ritz.py FCIDUMP OUTPUT_DIR")

    fcidump_path = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    nthreads = int(os.environ.get("OMP_NUM_THREADS", "64"))
    max_cycle = int(os.environ.get("FCI_MAX_CYCLE", "110"))
    max_space = int(os.environ.get("FCI_MAX_SPACE", "14"))
    max_memory_mb = int(os.environ.get("PYSCF_MAX_MEMORY", "450000"))
    lib.num_threads(nthreads)

    input_sha256 = sha256_file(fcidump_path)
    emit(
        "environment",
        host=platform.node(),
        python=sys.version.split()[0],
        numpy=np.__version__,
        pyscf=pyscf.__version__,
        threads=lib.num_threads(),
        max_cycle=max_cycle,
        max_space=max_space,
        max_memory_mb=max_memory_mb,
        fcidump=str(fcidump_path),
        fcidump_sha256=input_sha256,
        expected_fcidump_sha256=EXPECTED_FCIDUMP_SHA256,
    )
    if input_sha256 != EXPECTED_FCIDUMP_SHA256:
        raise RuntimeError("FCIDUMP SHA-256 does not match the frozen scientific input")

    data = fcidump.read(str(fcidump_path), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    ecore = float(data["ECORE"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    h1 = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)

    from math import comb

    nstr_alpha = comb(norb, nelec[0])
    nstr_beta = comb(norb, nelec[1])
    emit(
        "problem",
        norb=norb,
        nelec_total=nelec_total,
        ms2=ms2,
        nelec_alpha=nelec[0],
        nelec_beta=nelec[1],
        ecore=ecore,
        nstr_alpha=nstr_alpha,
        nstr_beta=nstr_beta,
        dimension=nstr_alpha * nstr_beta,
        vector_gib=nstr_alpha * nstr_beta * 8 / 2**30,
    )

    solver = fci.direct_spin1.FCI()
    solver.max_memory = max_memory_mb
    solver.max_cycle = max_cycle
    solver.max_space = max_space
    solver.conv_tol = 1e-10
    solver.verbose = 5
    solver.davidson_only = True
    solver.pspace_size = 800

    emit("kernel_start")
    kernel_started = time.time()
    energy_kernel, civec = solver.kernel(h1, eri, norb, nelec, ecore=ecore)
    kernel_wall_s = time.time() - kernel_started
    civec = np.asarray(civec, dtype=np.float64).reshape(nstr_alpha, nstr_beta)
    emit(
        "kernel_end",
        energy_kernel_eh=float(energy_kernel),
        solver_converged=bool(solver.converged),
        kernel_wall_s=kernel_wall_s,
        public_hci_bound_eh=PUBLIC_HCI_BOUND_EH,
        margin_eh=PUBLIC_HCI_BOUND_EH - float(energy_kernel),
    )

    vector_path = output_dir / "fci_ritz_vector.npy"
    np.save(vector_path, civec, allow_pickle=False)
    vector_sha256 = sha256_file(vector_path)
    emit(
        "vector_saved",
        path=str(vector_path),
        bytes=vector_path.stat().st_size,
        sha256=vector_sha256,
        norm2=float(np.vdot(civec.ravel(), civec.ravel()).real),
        max_abs_coefficient=float(np.max(np.abs(civec))),
    )

    # Fresh H|c> contraction after the eigensolver.  This recomputes both the
    # Rayleigh quotient and residual from the saved candidate, rather than
    # trusting the scalar printed by Davidson.
    emit("independent_contraction_start")
    contraction_started = time.time()
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, norb, nelec, 0.5)
    hc = np.asarray(
        fci.direct_spin1.contract_2e(h2e, civec, norb, nelec), dtype=np.float64
    ).reshape(civec.shape)
    norm2 = float(np.vdot(civec.ravel(), civec.ravel()).real)
    electronic_rayleigh = float(np.vdot(civec.ravel(), hc.ravel()).real / norm2)
    rayleigh = electronic_rayleigh + ecore
    residual = hc - electronic_rayleigh * civec
    residual_norm = float(np.linalg.norm(residual.ravel()) / np.sqrt(norm2))
    contraction_wall_s = time.time() - contraction_started

    result = {
        "schema_version": "2fe2s-fci-ritz-result/v2",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "pyscf_version": pyscf.__version__,
        "threads": lib.num_threads(),
        "max_cycle": max_cycle,
        "max_space": max_space,
        "conv_tol": solver.conv_tol,
        "solver_converged": bool(solver.converged),
        "fcidump_path": str(fcidump_path),
        "fcidump_sha256": input_sha256,
        "norb": norb,
        "nelec": list(nelec),
        "ms2": ms2,
        "dimension": nstr_alpha * nstr_beta,
        "ecore_eh": ecore,
        "energy_kernel_eh": float(energy_kernel),
        "energy_rayleigh_eh": rayleigh,
        "kernel_vs_rayleigh_abs_eh": abs(float(energy_kernel) - rayleigh),
        "residual_norm_eh": residual_norm,
        "norm2": norm2,
        "public_hci_bound_eh": PUBLIC_HCI_BOUND_EH,
        "improvement_vs_public_hci_eh": PUBLIC_HCI_BOUND_EH - rayleigh,
        "improvement_vs_public_hci_millihartree": 1000 * (PUBLIC_HCI_BOUND_EH - rayleigh),
        "vector_path": str(vector_path),
        "vector_bytes": vector_path.stat().st_size,
        "vector_sha256": vector_sha256,
        "kernel_wall_s": kernel_wall_s,
        "contraction_wall_s": contraction_wall_s,
        "total_wall_s": time.time() - START,
    }
    result_path = output_dir / "fci_ritz_result.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    emit("independent_contraction_end", **result)
    emit("complete", result_path=str(result_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
