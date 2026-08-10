#!/usr/bin/env python3
"""Build manuscript result macros from frozen, independently gated outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path


EXPECTED_GRID = [
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096,
    8192,
    16384,
    32768,
    65536,
    131072,
    262144,
    524288,
    1048576,
    2097152,
    4194304,
]


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def row_at(rows: list[dict[str, str]], k: int) -> dict[str, str]:
    matches = [row for row in rows if int(row["k"]) == k]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one K={k} row")
    return matches[0]


def latex_scientific(value: float, digits: int = 2) -> str:
    if value == 0.0:
        return "0"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10.0**exponent)
    return f"{mantissa:.{digits}f}\\times10^{{{exponent}}}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    verification = json.loads(
        (args.data_dir / "publication_state_verification.json").read_text(
            encoding="utf-8"
        )
    )
    benchmark = load_csv(args.data_dir / "topk_tracker_basis.csv")
    natural = load_csv(args.data_dir / "topk_natural_basis.csv")
    natural_summary = json.loads(
        (args.data_dir / "compressibility_natural_basis.json").read_text(
            encoding="utf-8"
        )
    )
    natural_rediagonalization = json.loads(
        (args.data_dir / "natural_top512_rediagonalization.json").read_text(
            encoding="utf-8"
        )
    )
    decomposition = load_csv(
        args.data_dir / "hamiltonian_decomposition_tracker_basis.csv"
    )
    if [int(row["k"]) for row in benchmark] != EXPECTED_GRID:
        raise SystemExit("benchmark top-K grid mismatch")
    if [int(row["k"]) for row in natural] != EXPECTED_GRID:
        raise SystemExit("natural top-K grid mismatch")
    if [int(row["k"]) for row in decomposition] != EXPECTED_GRID:
        raise SystemExit("decomposition grid mismatch")
    if not all(bool(value) for value in natural_summary["gates"].values()):
        raise SystemExit("natural-basis analysis contains a failed gate")
    if not all(
        bool(value)
        for value in natural_rediagonalization["gates"].values()
    ):
        raise SystemExit("natural-basis top-512 rediagonalization failed a gate")

    energy = float(verification["energy_rayleigh_eh"])
    residual = float(verification["residual_norm_eh"])
    if abs(energy + 116.60560912042631) > 2e-13:
        raise SystemExit("publication energy gate failed")
    if residual > 1e-6:
        raise SystemExit("publication residual gate failed")
    if abs(float(natural_summary["full_state_energy_eh"]) - energy) > 2e-9:
        raise SystemExit("natural-basis full-state energy gate failed")

    natural_trace = float(natural_summary["natural_occupation_trace"])
    natural_trace_expected = float(
        natural_summary["natural_occupation_trace_expected_from_norm"]
    )
    natural_trace_error = float(
        natural_summary["natural_occupation_trace_abs_error"]
    )
    natural_rdm1_symmetry = float(
        natural_summary["natural_rdm1_symmetry_max_abs"]
    )
    natural_rdm1_offdiagonal = float(
        natural_summary["natural_basis_rdm1_offdiagonal_max_abs"]
    )
    natural_rdm1_reference_difference = float(
        natural_summary["natural_basis_rdm1_reference_max_abs_difference"]
    )
    natural_occupations_reference_difference = float(
        natural_summary["reference_occupations_max_abs_difference"]
    )

    natural_fixed_512 = float(
        natural_rediagonalization["fixed_coefficient_energy_reference_eh"]
    )
    natural_relaxed_512 = float(
        natural_rediagonalization["same_support_optimized_energy_pyscf_eh"]
    )
    natural_relaxed_512_reference = float(
        natural_rediagonalization["same_support_optimized_energy_reference_eh"]
    )
    natural_relaxed_512_gain = float(
        natural_rediagonalization[
            "energy_gain_from_coefficient_reoptimization_eh"
        ]
    )
    natural_projected_residual = float(
        natural_rediagonalization["projected_residual_pyscf_eh"]
    )
    natural_full_residual = float(
        natural_rediagonalization["full_space_residual_pyscf_eh"]
    )
    if abs(natural_fixed_512 + 116.43662614444148) > 1e-10:
        raise SystemExit("natural-basis inherited top-512 energy gate failed")
    if natural_relaxed_512 > natural_fixed_512 + 1e-12:
        raise SystemExit("natural-basis top-512 relaxation is not variational")
    if abs(natural_relaxed_512 - natural_relaxed_512_reference) > 1e-10:
        raise SystemExit("natural-basis top-512 independent-energy gate failed")

    b_last = row_at(benchmark, 4194304)
    n_last = row_at(natural, 4194304)
    d_512 = row_at(decomposition, 512)
    weight = float(d_512["retained_norm2"])
    cross = float(d_512["retained_tail_cross_eh"])
    tail = float(d_512["tail_h_tail_eh"])
    decomposition_energy = float(d_512["full_energy_eh"])
    coupling_contribution = -cross / weight
    tail_contribution = -(
        tail - (1.0 - weight) * decomposition_energy
    ) / weight
    projected_error = float(d_512["projected_energy_error_eh"])
    if abs(coupling_contribution + tail_contribution - projected_error) > 2e-12:
        raise SystemExit("K=512 shifted decomposition gate failed")

    residual_last = float(n_last["residual_norm_eh"])
    residual_exponent = -2
    residual_mantissa = residual_last / (10.0**residual_exponent)
    final_residual_exponent = int(math.floor(math.log10(residual)))
    final_residual_mantissa = residual / (10.0**final_residual_exponent)
    text = (
        "% Frozen numerical values generated by analysis/build_results_tex.py.\n"
        f"\\newcommand{{\\FinalEnergy}}{{{energy:.14f}}}\n"
        f"\\newcommand{{\\FinalResidual}}{{{final_residual_mantissa:.4f}"
        f"\\times10^{{{final_residual_exponent}}}}}\n"
        f"\\newcommand{{\\BenchmarkKFiveTwelveCoupling}}{{{coupling_contribution:.7f}}}\n"
        f"\\newcommand{{\\BenchmarkKFiveTwelveTail}}{{{tail_contribution:.7f}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveRelaxedEnergy}}{{{natural_relaxed_512:.12f}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveRelaxedErrorMilli}}{{{1000.0 * (natural_relaxed_512 - energy):.3f}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveRelaxationGainMilli}}{{{1000.0 * natural_relaxed_512_gain:.3f}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveProjectedResidual}}{{{latex_scientific(natural_projected_residual)}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveFullResidual}}{{{natural_full_residual:.6f}}}\n"
        f"\\newcommand{{\\NaturalKFiveTwelveIndependentDifference}}{{{latex_scientific(abs(natural_relaxed_512 - natural_relaxed_512_reference))}}}\n"
        f"\\newcommand{{\\NaturalReplayTrace}}{{{natural_trace:.15f}}}\n"
        f"\\newcommand{{\\NaturalReplayExpectedTrace}}{{{natural_trace_expected:.15f}}}\n"
        f"\\newcommand{{\\NaturalReplayTraceError}}{{{latex_scientific(natural_trace_error)}}}\n"
        f"\\newcommand{{\\NaturalReplayRDMSymmetry}}{{{latex_scientific(natural_rdm1_symmetry)}}}\n"
        f"\\newcommand{{\\NaturalReplayRDMOffDiagonal}}{{{latex_scientific(natural_rdm1_offdiagonal)}}}\n"
        f"\\newcommand{{\\NaturalReplayRDMReferenceDifference}}{{{latex_scientific(natural_rdm1_reference_difference)}}}\n"
        f"\\newcommand{{\\NaturalReplayOccupationsDifference}}{{{latex_scientific(natural_occupations_reference_difference)}}}\n"
        "\\newcommand{\\LateBasisSentence}{At \\(K=4{,}194{,}304\\), "
        "the natural-orbital projection increases \\(W_K\\) from "
        f"{float(b_last['cumulative_norm2']):.6f} to "
        f"{float(n_last['cumulative_norm2']):.6f} and reduces "
        "\\(\\Delta E_K\\) from "
        f"{1000.0 * float(b_last['energy_error_vs_full_eh']):.3f} to "
        f"\\({1000.0 * float(n_last['energy_error_vs_full_eh']):.4f}"
        "\\,\\mathrm{m}\\Eh\\); its residual nevertheless remains "
        f"\\({residual_mantissa:.2f}\\times10^{{{residual_exponent}}}"
        "\\,\\Eh\\).}\n"
    )
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, args.output)
    print(text, end="")


if __name__ == "__main__":
    main()
