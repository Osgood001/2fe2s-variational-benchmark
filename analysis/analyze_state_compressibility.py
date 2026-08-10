#!/usr/bin/env python3
"""Basis-resolved determinant compressibility of a frozen CASCI Ritz state.

The full vector is validated by SHA-256, contracted directly with the FCIDUMP
Hamiltonian, and then projected onto the determinants with the largest squared
coefficients.  In addition to norm, energy, and residual curves, the program
decomposes the full-state energy for p=P_K c and q=(1-P_K)c as

    <p|H|p> + 2 Re <p|H|q> + <q|H|q>.

Only compact CSV/JSON results are written; the full vector remains an input
artifact.  This program is intended for a hash-locked Bohrium job.
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
from math import comb
from pathlib import Path

import numpy as np
from pyscf import __version__ as pyscf_version
from pyscf import ao2mo, fci, lib
from pyscf.tools import fcidump


DEFAULT_K_VALUES = (
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
)
RANK_CHUNK = 8_388_608
MOMENT_CHUNK = 8_388_608


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("fcidump", type=Path)
    parser.add_argument("vector", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--basis-label", required=True)
    parser.add_argument("--expected-fcidump-sha256", required=True)
    parser.add_argument("--expected-vector-sha256")
    parser.add_argument("--reference-energy-eh", type=float)
    parser.add_argument("--reference-energy-tolerance-eh", type=float, default=2e-10)
    parser.add_argument("--rhf-energy-eh", type=float)
    parser.add_argument("--max-full-residual-eh", type=float)
    parser.add_argument(
        "--occupation-trace-abs-tolerance", type=float, default=1e-8
    )
    parser.add_argument("--reference-occupations-csv", type=Path)
    parser.add_argument("--expected-reference-occupations-sha256")
    parser.add_argument(
        "--reference-occupations-max-abs-tolerance", type=float, default=2e-8
    )
    parser.add_argument(
        "--natural-basis-rdm1-max-abs-tolerance", type=float, default=2e-8
    )
    parser.add_argument(
        "--write-natural-orbital-fcidump",
        action="store_true",
        help="write the natural-orbital rotation and transformed Hamiltonian",
    )
    parser.add_argument(
        "--k-values",
        default=",".join(str(value) for value in DEFAULT_K_VALUES),
        help="comma-separated determinant counts",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def global_top_k(flat: np.ndarray, kmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Return exact global top-k indices and squared amplitudes, descending."""
    kept_indices = np.empty(0, dtype=np.int64)
    kept_weights = np.empty(0, dtype=np.float64)
    for start in range(0, flat.size, RANK_CHUNK):
        stop = min(start + RANK_CHUNK, flat.size)
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
        print(f"RANK progress={stop}/{flat.size} retained={kept_weights.size}", flush=True)
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


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_reference_occupations(path: Path, norb: int) -> np.ndarray:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != norb:
        raise SystemExit(
            f"reference-occupation row count mismatch: {len(rows)} != {norb}"
        )
    expected_indices = list(range(1, norb + 1))
    actual_indices = [int(row["natural_orbital"]) for row in rows]
    if actual_indices != expected_indices:
        raise SystemExit("reference occupations are not in canonical order")
    return np.asarray([float(row["occupation"]) for row in rows])


def main() -> None:
    args = parse_args()
    started = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "64")))
    k_values = tuple(int(token) for token in args.k_values.split(","))
    if sorted(set(k_values)) != list(k_values) or min(k_values) < 1:
        raise SystemExit("k-values must be positive, unique, and ascending")

    actual_reference_occupations_hash: str | None = None
    if args.reference_occupations_csv is not None:
        if args.expected_reference_occupations_sha256 is None:
            raise SystemExit(
                "--expected-reference-occupations-sha256 is required with "
                "--reference-occupations-csv"
            )
        actual_reference_occupations_hash = sha256_file(
            args.reference_occupations_csv
        )
        if (
            actual_reference_occupations_hash
            != args.expected_reference_occupations_sha256
        ):
            raise SystemExit(
                "reference-occupation CSV hash mismatch: "
                f"{actual_reference_occupations_hash}"
            )
    elif args.expected_reference_occupations_sha256 is not None:
        raise SystemExit(
            "--reference-occupations-csv is required with "
            "--expected-reference-occupations-sha256"
        )

    actual_fcidump_hash = sha256_file(args.fcidump)
    actual_vector_hash = sha256_file(args.vector)
    if actual_fcidump_hash != args.expected_fcidump_sha256:
        raise SystemExit(f"FCIDUMP hash mismatch: {actual_fcidump_hash}")
    if args.expected_vector_sha256 and actual_vector_hash != args.expected_vector_sha256:
        raise SystemExit(f"vector hash mismatch: {actual_vector_hash}")

    data = fcidump.read(str(args.fcidump), verbose=False)
    norb = int(data["NORB"])
    nelec_total = int(data["NELEC"])
    ms2 = int(data["MS2"])
    nelec = ((nelec_total + ms2) // 2, (nelec_total - ms2) // 2)
    ecore = float(data["ECORE"])
    h1e = np.asarray(data["H1"], dtype=np.float64)
    eri = ao2mo.restore(8, np.asarray(data["H2"], dtype=np.float64), norb)
    shape = (comb(norb, nelec[0]), comb(norb, nelec[1]))
    reference_occupations: np.ndarray | None = None
    if args.reference_occupations_csv is not None:
        reference_occupations = load_reference_occupations(
            args.reference_occupations_csv, norb
        )
    vector = np.load(args.vector, allow_pickle=False)
    if vector.shape != shape or vector.dtype != np.float64:
        raise SystemExit(
            f"unexpected vector: shape={vector.shape}, dtype={vector.dtype}; "
            f"expected {shape}, float64"
        )
    flat = vector.reshape(-1)
    if max(k_values) > flat.size:
        raise SystemExit("largest k exceeds Hilbert-space dimension")

    moments = coefficient_moments(flat)
    if abs(moments["norm2"] - 1.0) > 2e-12:
        raise SystemExit(f"vector norm gate failed: {moments['norm2']}")

    h2e = fci.direct_spin1.absorb_h1e(h1e, eri, norb, nelec, 0.5)
    full_contracted = np.asarray(
        fci.direct_spin1.contract_2e(h2e, vector, norb, nelec), dtype=np.float64
    ).reshape(shape)
    full_electronic_energy = float(np.dot(flat, full_contracted.reshape(-1)))
    full_energy = full_electronic_energy + ecore
    full_residual = float(
        np.linalg.norm((full_contracted - full_electronic_energy * vector).reshape(-1))
    )
    del full_contracted
    gc.collect()
    if (
        args.reference_energy_eh is not None
        and abs(full_energy - args.reference_energy_eh)
        > args.reference_energy_tolerance_eh
    ):
        raise SystemExit(
            f"full-state energy gate failed: {full_energy} vs {args.reference_energy_eh}"
        )

    top_indices, top_weights = global_top_k(flat, max(k_values))
    cumulative_weights = np.cumsum(top_weights, dtype=np.float64)
    rows: list[dict[str, object]] = []
    decomposition_rows: list[dict[str, object]] = []
    for k in k_values:
        k_started = time.time()
        retained_norm2 = float(cumulative_weights[k - 1])
        root_weight = float(np.sqrt(retained_norm2))
        projected = np.zeros_like(vector)
        projected_flat = projected.reshape(-1)
        selected = top_indices[:k]
        projected_flat[selected] = flat[selected]
        projected /= root_weight
        overlap = float(np.dot(flat, projected_flat))
        contracted = np.asarray(
            fci.direct_spin1.contract_2e(h2e, projected, norb, nelec),
            dtype=np.float64,
        ).reshape(shape)
        contracted_flat = contracted.reshape(-1)
        projected_electronic_energy = float(np.dot(projected_flat, contracted_flat))
        projected_energy = projected_electronic_energy + ecore
        c_h_phi = float(np.dot(flat, contracted_flat) + ecore * overlap)
        residual = float(
            np.linalg.norm(
                (contracted - projected_electronic_energy * projected).reshape(-1)
            )
        )

        # p=P_K c=sqrt(W_K) phi_K and q=c-p.  All terms below include ECORE.
        retained_term = retained_norm2 * projected_energy
        p_h_c = root_weight * c_h_phi
        retained_tail_half = p_h_c - retained_term
        cross_term = 2.0 * retained_tail_half
        tail_term = full_energy - retained_term - cross_term
        decomposition_closure = retained_term + cross_term + tail_term - full_energy

        # Equivalent decomposition of E_K-E_full in the raw H convention.
        # The three individual terms below depend on the arbitrary additive
        # constant in H, so they are retained only as an exact bookkeeping
        # cross-check.
        normalization_term = (1.0 - retained_norm2) * projected_energy
        removed_tail_term = -tail_term
        removed_cross_term = -cross_term
        energy_error = projected_energy - full_energy
        error_closure = (
            normalization_term + removed_tail_term + removed_cross_term - energy_error
        )

        # Gauge-invariant decomposition with A = H-E_full I.  Because
        # <c|A|c>=0, W_K(E_K-E_full) = -<q|A|q>-2<p|A|q>.  Dividing by W_K
        # separates the normalized projection error into a tail-quadratic
        # contribution and a retained--tail coupling contribution without any
        # dependence on the chosen Hamiltonian zero.
        retained_shifted_term = retained_term - retained_norm2 * full_energy
        tail_norm2 = 1.0 - retained_norm2
        tail_shifted_term = tail_term - tail_norm2 * full_energy
        shifted_closure = retained_shifted_term + cross_term + tail_shifted_term
        tail_quadratic_error_contribution = -tail_shifted_term / retained_norm2
        retained_tail_coupling_error_contribution = -cross_term / retained_norm2
        shifted_error_closure = (
            tail_quadratic_error_contribution
            + retained_tail_coupling_error_contribution
            - energy_error
        )
        recovery = None
        if args.rhf_energy_eh is not None:
            recovery = float(
                (args.rhf_energy_eh - projected_energy)
                / (args.rhf_energy_eh - full_energy)
            )

        row = {
            "basis": args.basis_label,
            "k": k,
            "fraction_of_full_space": k / flat.size,
            "cumulative_norm2": retained_norm2,
            "projected_overlap": overlap,
            "rayleigh_energy_eh": projected_energy,
            "energy_error_vs_full_eh": energy_error,
            "energy_error_vs_full_millieh": 1e3 * energy_error,
            "residual_norm_eh": residual,
            "correlation_recovery_fraction": recovery,
            "wall_s": time.time() - k_started,
            "state": "amplitude-ranked projection",
        }
        decomposition = {
            "basis": args.basis_label,
            "k": k,
            "retained_norm2": retained_norm2,
            "retained_h_retained_eh": retained_term,
            "retained_tail_cross_eh": cross_term,
            "tail_h_tail_eh": tail_term,
            "full_energy_eh": full_energy,
            "decomposition_closure_eh": decomposition_closure,
            "normalization_contribution_to_error_eh": normalization_term,
            "removed_tail_contribution_to_error_eh": removed_tail_term,
            "removed_cross_contribution_to_error_eh": removed_cross_term,
            "projected_energy_error_eh": energy_error,
            "energy_error_closure_eh": error_closure,
            "retained_shifted_h_retained_eh": retained_shifted_term,
            "tail_norm2": tail_norm2,
            "tail_shifted_h_tail_eh": tail_shifted_term,
            "shifted_decomposition_closure_eh": shifted_closure,
            "tail_quadratic_error_contribution_eh": tail_quadratic_error_contribution,
            "retained_tail_coupling_error_contribution_eh": retained_tail_coupling_error_contribution,
            "shifted_energy_error_closure_eh": shifted_error_closure,
        }
        rows.append(row)
        decomposition_rows.append(decomposition)
        print("TOPK_RESULT " + json.dumps(row, sort_keys=True), flush=True)
        print("DECOMPOSITION " + json.dumps(decomposition, sort_keys=True), flush=True)
        del contracted, contracted_flat, projected, projected_flat, selected
        gc.collect()

    rdm1 = fci.direct_spin1.make_rdm1(vector, norb, nelec)
    rdm1_symmetry_max_abs = float(np.max(np.abs(rdm1 - rdm1.T)))
    symmetrized_rdm1 = 0.5 * (rdm1 + rdm1.T)
    occupations, natural_orbitals = np.linalg.eigh(symmetrized_rdm1)
    order = np.argsort(occupations)[::-1]
    occupations = occupations[order]
    natural_orbitals = natural_orbitals[:, order]
    # Canonicalize the otherwise arbitrary sign of every real eigenvector.
    for column in range(norb):
        pivot = int(np.argmax(np.abs(natural_orbitals[:, column])))
        if natural_orbitals[pivot, column] < 0.0:
            natural_orbitals[:, column] *= -1.0
    occupation_trace = float(np.sum(occupations))
    expected_occupation_trace = float(nelec_total * moments["norm2"])
    occupation_trace_abs_error = abs(
        occupation_trace - expected_occupation_trace
    )
    occupation_trace_relative_error = (
        occupation_trace_abs_error / abs(expected_occupation_trace)
    )
    reference_occupations_max_abs_difference: float | None = None
    natural_basis_rdm1_offdiagonal_max_abs = float(
        np.max(
            np.abs(
                symmetrized_rdm1
                - np.diag(np.diag(symmetrized_rdm1))
            )
        )
    )
    natural_basis_rdm1_reference_max_abs_difference: float | None = None
    if reference_occupations is not None:
        reference_occupations_max_abs_difference = float(
            np.max(np.abs(occupations - reference_occupations))
        )
        natural_basis_rdm1_reference_max_abs_difference = float(
            np.max(
                np.abs(
                    symmetrized_rdm1 - np.diag(reference_occupations)
                )
            )
        )
    occupation_rows = [
        {
            "basis": args.basis_label,
            "natural_orbital": orbital,
            "occupation": float(occupation),
            "distance_to_idempotent": float(
                min(abs(occupation), abs(2.0 - occupation))
            ),
        }
        for orbital, occupation in enumerate(occupations, start=1)
    ]

    derived_natural_basis: dict[str, object] | None = None
    if args.write_natural_orbital_fcidump:
        rotation_path = args.output_dir / "natural_orbital_rotation.csv"
        with rotation_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["source_orbital", "natural_orbital", "coefficient"])
            for source_orbital in range(norb):
                for natural_orbital in range(norb):
                    writer.writerow(
                        [
                            source_orbital + 1,
                            natural_orbital + 1,
                            format(
                                float(natural_orbitals[source_orbital, natural_orbital]),
                                ".17g",
                            ),
                        ]
                    )
        natural_h1e = natural_orbitals.T @ h1e @ natural_orbitals
        natural_eri = ao2mo.incore.full(eri, natural_orbitals, compact=True)
        natural_fcidump_path = args.output_dir / "natural_orbital.fcidump"
        fcidump.from_integrals(
            str(natural_fcidump_path),
            natural_h1e,
            natural_eri,
            norb,
            nelec_total,
            nuc=ecore,
            ms=ms2,
            orbsym=[1] * norb,
            tol=1e-15,
            float_format=" %.17e",
        )
        rotated_rdm1 = natural_orbitals.T @ rdm1 @ natural_orbitals
        derived_natural_basis = {
            "fcidump_filename": natural_fcidump_path.name,
            "fcidump_sha256": sha256_file(natural_fcidump_path),
            "fcidump_bytes": natural_fcidump_path.stat().st_size,
            "rotation_filename": rotation_path.name,
            "rotation_sha256": sha256_file(rotation_path),
            "rotation_bytes": rotation_path.stat().st_size,
            "rotation_orthogonality_max_abs": float(
                np.max(np.abs(natural_orbitals.T @ natural_orbitals - np.eye(norb)))
            ),
            "rotated_rdm1_offdiagonal_max_abs": float(
                np.max(np.abs(rotated_rdm1 - np.diag(np.diag(rotated_rdm1))))
            ),
            "source_fcidump_sha256": actual_fcidump_hash,
            "source_vector_sha256": actual_vector_hash,
        }

    write_csv(args.output_dir / "topk_compressibility.csv", rows)
    write_csv(args.output_dir / "hamiltonian_decomposition.csv", decomposition_rows)
    write_csv(args.output_dir / "natural_occupations.csv", occupation_rows)
    gates = {
        "fcidump_hash": actual_fcidump_hash == args.expected_fcidump_sha256,
        "vector_hash": not args.expected_vector_sha256
        or actual_vector_hash == args.expected_vector_sha256,
        "reference_occupations_hash": (
            args.reference_occupations_csv is None
            or actual_reference_occupations_hash
            == args.expected_reference_occupations_sha256
        ),
        "norm": abs(moments["norm2"] - 1.0) <= 2e-12,
        "reference_energy": args.reference_energy_eh is None
        or abs(full_energy - args.reference_energy_eh)
        <= args.reference_energy_tolerance_eh,
        "full_residual": args.max_full_residual_eh is None
        or full_residual <= args.max_full_residual_eh,
        "natural_occupation_trace": occupation_trace_abs_error
        <= args.occupation_trace_abs_tolerance,
        "natural_occupation_bounds": float(np.min(occupations)) >= -2e-8
        and float(np.max(occupations)) <= 2.0 + 2e-8,
        "natural_rdm1_symmetry": rdm1_symmetry_max_abs <= 2e-8,
        "natural_occupations_reference": (
            reference_occupations_max_abs_difference is None
            or reference_occupations_max_abs_difference
            <= args.reference_occupations_max_abs_tolerance
        ),
        "natural_basis_rdm1_offdiagonal": (
            reference_occupations is None
            or natural_basis_rdm1_offdiagonal_max_abs
            <= args.natural_basis_rdm1_max_abs_tolerance
        ),
        "natural_basis_rdm1_reference": (
            natural_basis_rdm1_reference_max_abs_difference is None
            or natural_basis_rdm1_reference_max_abs_difference
            <= args.natural_basis_rdm1_max_abs_tolerance
        ),
        "topk_variational": all(
            float(row["rayleigh_energy_eh"]) >= full_energy - 1e-10 for row in rows
        ),
        "topk_nested_norm": all(
            float(rows[index]["cumulative_norm2"])
            <= float(rows[index + 1]["cumulative_norm2"])
            for index in range(len(rows) - 1)
        ),
        "decomposition_closure": max(
            abs(float(row["decomposition_closure_eh"]))
            for row in decomposition_rows
        )
        <= 2e-12,
        "energy_error_closure": max(
            abs(float(row["energy_error_closure_eh"]))
            for row in decomposition_rows
        )
        <= 2e-12,
        "shifted_decomposition_closure": max(
            abs(float(row["shifted_decomposition_closure_eh"]))
            for row in decomposition_rows
        )
        <= 2e-12,
        "shifted_energy_error_closure": max(
            abs(float(row["shifted_energy_error_closure_eh"]))
            for row in decomposition_rows
        )
        <= 2e-12,
    }
    summary = {
        "schema_version": "2fe2s-basis-compressibility/v3",
        "basis": args.basis_label,
        "generated_at_unix": time.time(),
        "wall_s": time.time() - started,
        "host": platform.node(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pyscf_version": pyscf_version,
        "threads": lib.num_threads(),
        "fcidump_sha256": actual_fcidump_hash,
        "vector_sha256": actual_vector_hash,
        "vector_bytes": args.vector.stat().st_size,
        "vector_shape": list(vector.shape),
        "dimension": int(flat.size),
        "norb": norb,
        "nelec": list(nelec),
        "ms2": ms2,
        "ecore_eh": ecore,
        "full_state_energy_eh": full_energy,
        "full_state_residual_eh": full_residual,
        "coefficient_moments": moments,
        "top_k_max": max(k_values),
        "top_k_max_cumulative_norm2": float(cumulative_weights[-1]),
        "top_k_boundary_weight": float(top_weights[-1]),
        "natural_occupation_trace": occupation_trace,
        "natural_occupation_trace_expected_from_norm": expected_occupation_trace,
        "natural_occupation_trace_abs_error": occupation_trace_abs_error,
        "natural_occupation_trace_relative_error": occupation_trace_relative_error,
        "natural_occupation_trace_abs_tolerance": args.occupation_trace_abs_tolerance,
        "natural_rdm1_symmetry_max_abs": rdm1_symmetry_max_abs,
        "reference_occupations_csv": (
            str(args.reference_occupations_csv)
            if args.reference_occupations_csv
            else None
        ),
        "reference_occupations_sha256": actual_reference_occupations_hash,
        "expected_reference_occupations_sha256": (
            args.expected_reference_occupations_sha256
        ),
        "reference_occupations_max_abs_difference": reference_occupations_max_abs_difference,
        "reference_occupations_max_abs_tolerance": args.reference_occupations_max_abs_tolerance,
        "natural_basis_rdm1_offdiagonal_max_abs": (
            natural_basis_rdm1_offdiagonal_max_abs
        ),
        "natural_basis_rdm1_reference_max_abs_difference": (
            natural_basis_rdm1_reference_max_abs_difference
        ),
        "natural_basis_rdm1_max_abs_tolerance": (
            args.natural_basis_rdm1_max_abs_tolerance
        ),
        "natural_occupations_descending": [float(value) for value in occupations],
        "topk_rows": rows,
        "hamiltonian_decomposition_rows": decomposition_rows,
        "derived_natural_basis": derived_natural_basis,
        "gates": gates,
    }
    (args.output_dir / "compressibility_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("SUMMARY " + json.dumps(gates, sort_keys=True), flush=True)
    if not all(gates.values()):
        raise SystemExit("one or more analysis gates failed: " + json.dumps(gates))


if __name__ == "__main__":
    main()
