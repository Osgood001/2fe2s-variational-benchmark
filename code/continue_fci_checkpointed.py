#!/usr/bin/env python3
"""Checkpointed warm-start Davidson continuation for 2fe_2s_30e_20o.

The input and every reported candidate are explicit CI vectors.  Energies and
residuals are recomputed by direct Hamiltonian contraction; convergence is not
inferred from a printed scalar.  Periodic callback checkpoints make long runs
recoverable without changing PySCF's Davidson algorithm.
"""

from __future__ import annotations

import gc
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


EXPECTED_FCIDUMP_SHA256 = os.environ.get(
    "EXPECTED_FCIDUMP_SHA256",
    "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7",
)
PUBLIC_HCI_BOUND_EH = -116.605425
START = time.time()


def sha256_file(path: Path, block_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def atomic_npy(path: Path, array: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def emit(event: str, **values: object) -> None:
    payload = {"elapsed_s": round(time.time() - START, 3), "event": event, **values}
    print("SIMULATION_JSON " + json.dumps(payload, sort_keys=True), flush=True)


def append_jsonl(path: Path, value: object) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(value, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def diagnose(
    civec: np.ndarray,
    h2e: np.ndarray,
    norb: int,
    nelec: tuple[int, int],
    ecore: float,
) -> dict[str, float]:
    started = time.time()
    hc = np.asarray(
        fci.direct_spin1.contract_2e(h2e, civec, norb, nelec), dtype=np.float64
    ).reshape(civec.shape)
    norm2 = float(np.vdot(civec.ravel(), civec.ravel()).real)
    electronic_energy = float(np.vdot(civec.ravel(), hc.ravel()).real / norm2)
    residual_norm = float(
        np.linalg.norm((hc - electronic_energy * civec).ravel()) / np.sqrt(norm2)
    )
    del hc
    gc.collect()
    return {
        "energy_rayleigh_eh": electronic_energy + ecore,
        "electronic_energy_eh": electronic_energy,
        "norm2": norm2,
        "residual_norm_eh": residual_norm,
        "diagnostic_wall_s": time.time() - started,
    }


def scalar(value: object) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: continue_fci_checkpointed.py FCIDUMP INITIAL_VECTOR_NPY OUTPUT_DIR"
        )

    fcidump_path = Path(sys.argv[1]).resolve()
    initial_vector_path = Path(sys.argv[2]).resolve()
    output_dir = Path(sys.argv[3]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    threads = int(os.environ.get("OMP_NUM_THREADS", "64"))
    max_cycle = int(os.environ.get("FCI_MAX_CYCLE", "160"))
    max_space = int(os.environ.get("FCI_MAX_SPACE", "64"))
    max_memory_mb = int(os.environ.get("PYSCF_MAX_MEMORY", "450000"))
    checkpoint_every = int(os.environ.get("FCI_CHECKPOINT_EVERY", "8"))
    energy_tolerance = float(os.environ.get("FCI_ENERGY_TOL", "1e-13"))
    residual_tolerance = float(os.environ.get("FCI_RESIDUAL_TOL", "1e-8"))
    lib.num_threads(threads)

    fcidump_sha256 = sha256_file(fcidump_path)
    initial_vector_sha256 = sha256_file(initial_vector_path)
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
        checkpoint_every=checkpoint_every,
        energy_tolerance=energy_tolerance,
        residual_tolerance=residual_tolerance,
        fcidump_sha256=fcidump_sha256,
        initial_vector_sha256=initial_vector_sha256,
    )
    if fcidump_sha256 != EXPECTED_FCIDUMP_SHA256:
        raise RuntimeError("FCIDUMP SHA-256 mismatch")

    data = fcidump.read(str(fcidump_path), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    ecore = float(data["ECORE"])
    h1 = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, norb, nelec, 0.5)

    from math import comb

    shape = (comb(norb, nelec[0]), comb(norb, nelec[1]))
    civec = np.load(initial_vector_path, allow_pickle=False)
    if civec.shape != shape or civec.dtype != np.float64:
        raise RuntimeError(
            f"initial vector shape/dtype mismatch: {civec.shape} {civec.dtype} != {shape} float64"
        )

    initial_diagnostics = diagnose(civec, h2e, norb, nelec, ecore)
    initial_record = {
        "schema_version": "2fe2s-fci-continuation-initial/v1",
        "fcidump_sha256": fcidump_sha256,
        "initial_vector_path": str(initial_vector_path),
        "initial_vector_sha256": initial_vector_sha256,
        "vector_shape": list(shape),
        **initial_diagnostics,
    }
    atomic_json(output_dir / "initial_diagnostics.json", initial_record)
    emit("initial_diagnostics", **initial_record)

    trajectory_path = output_dir / "continuation_trajectory.jsonl"
    if trajectory_path.exists():
        raise RuntimeError(f"refusing to overwrite existing trajectory: {trajectory_path}")

    checkpoint_path = output_dir / "latest_checkpoint.npy"
    checkpoint_metadata_path = output_dir / "latest_checkpoint.json"
    callback_saves = 0

    def callback(envs: dict[str, object]) -> None:
        nonlocal callback_saves
        cycle = int(envs["icyc"])
        energy_electronic = scalar(envs["e"])
        residual_norm = scalar(envs["dx_norm"])
        delta_energy = scalar(envs["de"])
        record = {
            "elapsed_s": time.time() - START,
            "cycle": cycle,
            "space": int(envs["space"]),
            "energy_eh": energy_electronic + ecore,
            "delta_energy_eh": delta_energy,
            "residual_norm_eh": residual_norm,
            "fresh_start_next": bool(envs["fresh_start"]),
        }
        append_jsonl(trajectory_path, record)
        if checkpoint_every > 0 and (cycle + 1) % checkpoint_every == 0:
            vector = np.asarray(envs["x0"][0], dtype=np.float64).reshape(shape)
            save_started = time.time()
            atomic_npy(checkpoint_path, vector)
            vector_sha256 = sha256_file(checkpoint_path)
            callback_saves += 1
            checkpoint_record = {
                "schema_version": "2fe2s-fci-callback-checkpoint/v1",
                **record,
                "checkpoint_every": checkpoint_every,
                "checkpoint_saves": callback_saves,
                "vector_path": str(checkpoint_path),
                "vector_bytes": checkpoint_path.stat().st_size,
                "vector_sha256": vector_sha256,
                "save_and_hash_wall_s": time.time() - save_started,
            }
            atomic_json(checkpoint_metadata_path, checkpoint_record)
            emit("checkpoint_saved", **checkpoint_record)

    solver = fci.direct_spin1.FCI()
    solver.max_memory = max_memory_mb
    solver.max_cycle = max_cycle
    solver.max_space = max_space
    solver.conv_tol = energy_tolerance
    solver.conv_tol_residual = residual_tolerance
    # PySCF 2.14 consumes this optional attribute in kernel_ms1 but does not
    # declare it in FCISolver._keys, so register it to keep sanity-check logs clean.
    solver._keys = solver._keys.union({"conv_tol_residual"})
    solver.verbose = 5
    solver.davidson_only = True
    solver.pspace_size = 800

    emit("kernel_start")
    kernel_started = time.time()
    energy_kernel, final_vector = solver.kernel(
        h1, eri, norb, nelec, ci0=civec, ecore=ecore, callback=callback
    )
    kernel_wall_s = time.time() - kernel_started
    final_vector = np.asarray(final_vector, dtype=np.float64).reshape(shape)
    emit(
        "kernel_end",
        energy_kernel_eh=float(energy_kernel),
        solver_converged=bool(solver.converged),
        kernel_wall_s=kernel_wall_s,
    )

    final_vector_path = output_dir / "fci_ritz_vector.npy"
    atomic_npy(final_vector_path, final_vector)
    final_vector_sha256 = sha256_file(final_vector_path)
    emit(
        "final_vector_saved",
        vector_path=str(final_vector_path),
        vector_bytes=final_vector_path.stat().st_size,
        vector_sha256=final_vector_sha256,
    )

    final_diagnostics = diagnose(final_vector, h2e, norb, nelec, ecore)
    result = {
        "schema_version": "2fe2s-fci-checkpointed-continuation/v1",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "pyscf_version": pyscf.__version__,
        "threads": lib.num_threads(),
        "max_cycle": max_cycle,
        "max_space": max_space,
        "max_memory_mb": max_memory_mb,
        "checkpoint_every": checkpoint_every,
        "callback_saves": callback_saves,
        "conv_tol_energy": energy_tolerance,
        "conv_tol_residual": residual_tolerance,
        "solver_converged": bool(solver.converged),
        "fcidump_path": str(fcidump_path),
        "fcidump_sha256": fcidump_sha256,
        "initial_vector_path": str(initial_vector_path),
        "initial_vector_sha256": initial_vector_sha256,
        "initial_energy_rayleigh_eh": initial_diagnostics["energy_rayleigh_eh"],
        "initial_residual_norm_eh": initial_diagnostics["residual_norm_eh"],
        "norb": norb,
        "nelec": list(nelec),
        "ms2": ms2,
        "dimension": shape[0] * shape[1],
        "ecore_eh": ecore,
        "energy_kernel_eh": float(energy_kernel),
        **final_diagnostics,
        "kernel_vs_rayleigh_abs_eh": abs(
            float(energy_kernel) - final_diagnostics["energy_rayleigh_eh"]
        ),
        "improvement_vs_public_hci_eh": PUBLIC_HCI_BOUND_EH
        - final_diagnostics["energy_rayleigh_eh"],
        "energy_drop_from_initial_eh": initial_diagnostics["energy_rayleigh_eh"]
        - final_diagnostics["energy_rayleigh_eh"],
        "vector_path": str(final_vector_path),
        "vector_bytes": final_vector_path.stat().st_size,
        "vector_sha256": final_vector_sha256,
        "kernel_wall_s": kernel_wall_s,
        "total_wall_s": time.time() - START,
    }
    atomic_json(output_dir / "continuation_result.json", result)
    emit("complete", **result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
