#!/usr/bin/env python3
"""Measure norm and energy compressibility of a hash-frozen CASCI Ritz state.

This program is intended to run as a Bohrium job.  It never writes the full
coefficient vector to its output; only compact CSV/JSON summaries are emitted.
"""

from __future__ import annotations

import csv
import gc
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
from pyscf import __version__ as pyscf_version
from pyscf import ao2mo, fci, lib
from pyscf.tools import fcidump


FCIDUMP_SHA256 = "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7"
VECTOR_SHA256 = "45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945"
FULL_ENERGY_EH = -116.60560912042631
FULL_RESIDUAL_EH = 6.901557706123275e-7
RHF_ENERGY_EH = -116.20581611838442
HCI_ENERGY_EH = -116.605425
DMRG_ROUNDED_EH = -116.6056091
COMPETITION_K512_ENERGY_EH = -116.3704626758453

K_VALUES = [
    64,
    128,
    256,
    512,
    1_024,
    2_048,
    4_096,
    8_192,
    16_384,
    32_768,
    65_536,
    131_072,
    262_144,
    524_288,
    1_048_576,
    2_097_152,
    4_194_304,
]
TOP_K_MAX = max(K_VALUES)
RANK_CHUNK = 8_388_608
MOMENT_CHUNK = 8_388_608


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def global_top_k(flat: np.ndarray, kmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Return exact global top-k indices and squared amplitudes, descending."""
    kept_indices = np.empty(0, dtype=np.int64)
    kept_weights = np.empty(0, dtype=np.float64)
    n = flat.size
    for start in range(0, n, RANK_CHUNK):
        stop = min(start + RANK_CHUNK, n)
        block = np.asarray(flat[start:stop])
        weights = np.square(block, dtype=np.float64)
        local_k = min(kmax, weights.size)
        if local_k < weights.size:
            local_selection = np.argpartition(weights, weights.size - local_k)[-local_k:]
        else:
            local_selection = np.arange(weights.size, dtype=np.int64)
        local_weights = weights[local_selection]
        local_indices = local_selection.astype(np.int64, copy=False) + start

        if kept_weights.size:
            merged_weights = np.concatenate((kept_weights, local_weights))
            merged_indices = np.concatenate((kept_indices, local_indices))
        else:
            merged_weights = local_weights
            merged_indices = local_indices
        keep = min(kmax, merged_weights.size)
        if keep < merged_weights.size:
            selection = np.argpartition(merged_weights, merged_weights.size - keep)[-keep:]
            kept_weights = merged_weights[selection]
            kept_indices = merged_indices[selection]
        else:
            kept_weights = merged_weights
            kept_indices = merged_indices
        print(
            f"RANK progress={stop}/{n} retained={kept_weights.size}",
            flush=True,
        )
        del block, weights, local_selection, local_weights, local_indices
        if "merged_weights" in locals():
            del merged_weights, merged_indices
        gc.collect()

    order = np.argsort(kept_weights, kind="stable")[::-1]
    return kept_indices[order], kept_weights[order]


def coefficient_moments(flat: np.ndarray) -> dict[str, float]:
    norm2 = 0.0
    sum_p2 = 0.0
    entropy = 0.0
    max_abs = 0.0
    for start in range(0, flat.size, MOMENT_CHUNK):
        stop = min(start + MOMENT_CHUNK, flat.size)
        block = np.asarray(flat[start:stop])
        probabilities = np.square(block, dtype=np.float64)
        norm2 += float(np.sum(probabilities, dtype=np.float64))
        sum_p2 += float(np.dot(probabilities, probabilities))
        nonzero = probabilities[probabilities > 0.0]
        entropy -= float(np.dot(nonzero, np.log(nonzero)))
        max_abs = max(max_abs, float(np.max(np.abs(block))))
        del block, probabilities, nonzero
    return {
        "norm2": norm2,
        "inverse_participation_ratio": sum_p2,
        "effective_count_ipr": 1.0 / sum_p2,
        "shannon_entropy_nats": entropy,
        "effective_count_shannon": float(np.exp(entropy)),
        "max_abs_coefficient": max_abs,
    }


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: analyze_compressibility.py FCIDUMP VECTOR_NPY OUTPUT_DIR"
        )
    fcidump_path = Path(sys.argv[1])
    vector_path = Path(sys.argv[2])
    output_dir = Path(sys.argv[3])
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "64")))

    actual_fcidump_hash = sha256_file(fcidump_path)
    actual_vector_hash = sha256_file(vector_path)
    if actual_fcidump_hash != FCIDUMP_SHA256:
        raise SystemExit(f"FCIDUMP hash mismatch: {actual_fcidump_hash}")
    if actual_vector_hash != VECTOR_SHA256:
        raise SystemExit(f"vector hash mismatch: {actual_vector_hash}")

    data = fcidump.read(str(fcidump_path), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    nelec_tuple = (
        (nelec_total + ms2) // 2,
        (nelec_total - ms2) // 2,
    )
    ecore = float(data["ECORE"])
    h1e = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)
    if norb != 20 or nelec_tuple != (15, 15) or abs(ecore) > 1e-15:
        raise SystemExit(
            f"unexpected sector: norb={norb}, nelec={nelec_tuple}, ecore={ecore}"
        )
    vector = np.load(vector_path, allow_pickle=False)
    if vector.shape != (15_504, 15_504) or vector.dtype != np.float64:
        raise SystemExit(f"unexpected vector: shape={vector.shape}, dtype={vector.dtype}")
    flat = vector.reshape(-1)

    moments = coefficient_moments(flat)
    if abs(moments["norm2"] - 1.0) > 2e-12:
        raise SystemExit(f"vector norm gate failed: {moments['norm2']}")

    top_indices, top_weights = global_top_k(flat, TOP_K_MAX)
    cumulative_weights = np.cumsum(top_weights, dtype=np.float64)
    boundary_gap = float(top_weights[-1])
    print(
        f"TOPK exact k={TOP_K_MAX} cumulative_norm2={cumulative_weights[-1]:.16g} "
        f"boundary_weight={boundary_gap:.16g}",
        flush=True,
    )

    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, nelec_tuple, 0.5)
    rows: list[dict[str, float | int | str]] = []
    for k in K_VALUES:
        k_started = time.time()
        norm2_before = float(cumulative_weights[k - 1])
        projected = np.zeros_like(vector)
        projected_flat = projected.reshape(-1)
        selected = top_indices[:k]
        projected_flat[selected] = flat[selected]
        projected /= np.sqrt(norm2_before)
        overlap = float(np.dot(flat, projected_flat))
        contracted = fci.direct_spin1.contract_2e(
            h2e, projected, norb, nelec_tuple
        )
        electronic_energy = float(
            np.dot(projected_flat, np.asarray(contracted).reshape(-1))
        )
        energy = electronic_energy + float(ecore)
        contracted -= electronic_energy * projected
        residual = float(np.linalg.norm(contracted.reshape(-1)))
        recovery = float(
            (RHF_ENERGY_EH - energy) / (RHF_ENERGY_EH - FULL_ENERGY_EH)
        )
        row = {
            "k": k,
            "cumulative_norm2": norm2_before,
            "projected_overlap": overlap,
            "rayleigh_energy_eh": energy,
            "energy_error_vs_full_eh": energy - FULL_ENERGY_EH,
            "energy_error_vs_full_millieh": 1e3 * (energy - FULL_ENERGY_EH),
            "residual_norm_eh": residual,
            "correlation_recovery_fraction": recovery,
            "wall_s": time.time() - k_started,
            "state": "top-k projection",
        }
        rows.append(row)
        print("TOPK_RESULT " + json.dumps(row, sort_keys=True), flush=True)
        del contracted, projected, projected_flat, selected
        gc.collect()

    rdm1 = fci.direct_spin1.make_rdm1(vector, norb, nelec_tuple)
    occupations = np.linalg.eigvalsh(0.5 * (rdm1 + rdm1.T))[::-1]
    if abs(float(np.sum(occupations)) - 30.0) > 1e-8:
        raise SystemExit(f"natural-occupation trace failed: {np.sum(occupations)}")

    with (output_dir / "topk_compressibility.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with (output_dir / "natural_occupations.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["natural_orbital", "occupation", "distance_to_idempotent"]
        )
        writer.writeheader()
        for orbital, occupation in enumerate(occupations, start=1):
            writer.writerow(
                {
                    "natural_orbital": orbital,
                    "occupation": float(occupation),
                    "distance_to_idempotent": float(
                        min(abs(occupation), abs(2.0 - occupation))
                    ),
                }
            )

    summary = {
        "schema_version": "2fe2s-compressibility/v1",
        "generated_at_unix": time.time(),
        "wall_s": time.time() - started,
        "host": platform.node(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": pyscf_version,
        "threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        "fcidump_sha256": actual_fcidump_hash,
        "vector_sha256": actual_vector_hash,
        "dimension": int(flat.size),
        "norb": int(norb),
        "nelec": list(nelec_tuple),
        "ecore_eh": float(ecore),
        "full_state_energy_eh": FULL_ENERGY_EH,
        "full_state_residual_eh": FULL_RESIDUAL_EH,
        "rhf_energy_eh": RHF_ENERGY_EH,
        "official_hci_energy_eh": HCI_ENERGY_EH,
        "literature_dmrg_rounded_energy_eh": DMRG_ROUNDED_EH,
        "competition_optimized_k512_energy_eh": COMPETITION_K512_ENERGY_EH,
        "coefficient_moments": moments,
        "top_k_max": TOP_K_MAX,
        "top_k_max_cumulative_norm2": float(cumulative_weights[-1]),
        "top_k_boundary_weight": boundary_gap,
        "natural_occupations_descending": [float(x) for x in occupations],
        "topk_rows": rows,
        "gates": {
            "fcidump_hash": actual_fcidump_hash == FCIDUMP_SHA256,
            "vector_hash": actual_vector_hash == VECTOR_SHA256,
            "norm": abs(moments["norm2"] - 1.0) <= 2e-12,
            "natural_occupation_trace": abs(float(np.sum(occupations)) - 30.0) <= 1e-8,
            "topk_variational": all(row["rayleigh_energy_eh"] >= FULL_ENERGY_EH - 1e-10 for row in rows),
            "topk_nested_norm": all(
                rows[i]["cumulative_norm2"] <= rows[i + 1]["cumulative_norm2"]
                for i in range(len(rows) - 1)
            ),
        },
    }
    if not all(summary["gates"].values()):
        raise SystemExit("one or more compressibility gates failed")
    (output_dir / "compressibility_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("SUMMARY " + json.dumps(summary["gates"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
