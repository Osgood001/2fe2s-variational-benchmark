#!/usr/bin/env python3
"""Build the two manuscript figures from frozen compact data only."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
GREY = "#666666"
LIGHT_GREY = "#D9D9D9"
BLACK = "#111111"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Nimbus Roman", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.0,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "axes.linewidth": 0.7,
        "lines.linewidth": 1.25,
        "lines.markersize": 3.8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)


def read_csv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def gate_topk_csv(rows, summary, label):
    summary_rows = summary["topk_rows"]
    if len(rows) != len(summary_rows):
        raise SystemExit(f"{label} CSV/summary row-count mismatch")
    numeric_fields = (
        "cumulative_norm2",
        "rayleigh_energy_eh",
        "energy_error_vs_full_eh",
        "residual_norm_eh",
    )
    for csv_row, summary_row in zip(rows, summary_rows):
        if int(csv_row["k"]) != int(summary_row["k"]):
            raise SystemExit(f"{label} CSV/summary K mismatch")
        for field in numeric_fields:
            if abs(float(csv_row[field]) - float(summary_row[field])) > 5e-13:
                raise SystemExit(
                    f"{label} CSV/summary mismatch for K={csv_row['k']} "
                    f"field={field}"
                )


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3.0, pad=2.0)


def panel_label(ax, label, x=-0.18):
    ax.text(
        x,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=8.5,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


def save(fig, stem):
    fig.savefig(OUT / (stem + ".pdf"))
    fig.savefig(OUT / (stem + ".svg"))
    fig.savefig(OUT / (stem + ".png"), dpi=450)
    plt.close(fig)


def figure_one():
    rows = read_csv(DATA / "continuation.csv")
    cycle = np.array([float(row["cycle"]) for row in rows])
    energy = np.array([float(row["energy_eh"]) for row in rows])
    residual = np.array([float(row["residual_norm_eh"]) for row in rows])
    final_energy = -116.60560912042631
    hci = -116.605425
    dmrg = -116.6056091

    fig, axes = plt.subplots(1, 3, figsize=(7.05, 2.18))
    fig.subplots_adjust(left=0.072, right=0.992, bottom=0.25, top=0.91, wspace=0.39)

    ax = axes[0]
    methods = ["HCI\n(reported)", "DMRG\n(rounded)", "Explicit Ritz\n(this work)"]
    y = np.array([2.0, 1.0, 0.0])
    x = 1e6 * (np.array([hci, dmrg, final_energy]) - hci)
    ax.hlines(y[1], x[1] - 0.05, x[1] + 0.05, color=ORANGE, lw=2.0)
    ax.plot(x[1], y[1], "s", ms=4.2, mfc="white", mec=ORANGE, mew=1.0)
    ax.plot(x[0], y[0], "D", color=BLACK, ms=4.2)
    ax.plot(x[2], y[2], "o", color=BLUE, ms=4.7)
    ax.axvline(0.0, color=GREY, lw=0.7, ls=(0, (2, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels(methods)
    ax.set_xlim(-205.0, 12.0)
    ax.set_ylim(-0.55, 2.55)
    ax.set_xlabel(r"$(E-E_{\rm HCI})\;(\mu E_{\rm h})$")
    ax.text(
        x[2] + 5.0,
        y[2],
        r"$\approx-184$",
        va="center",
        color=BLUE,
        fontsize=8.0,
    )
    ax.text(
        x[1] + 5.0,
        y[1],
        "published rounding interval",
        va="center",
        color=ORANGE,
        fontsize=7.4,
    )
    clean_axis(ax)
    panel_label(ax, "(a)", x=-0.27)

    ax = axes[1]
    callback = (cycle >= 0) & (cycle < np.max(cycle))
    gap_micro = 1e6 * (energy[callback] - final_energy)
    positive = gap_micro > 0.0
    ax.semilogy(
        cycle[callback][positive],
        gap_micro[positive],
        color=BLUE,
        marker="o",
        markevery=16,
        mfc="white",
        mec=BLUE,
        mew=0.8,
    )
    ax.set_xlabel("Davidson update")
    ax.set_ylabel(r"$E_m-E_{\rm Ritz}\;(\mu E_{\rm h})$")
    ax.set_xlim(0, 160)
    ax.set_ylim(1e-7, 1e2)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.yaxis.set_minor_formatter(NullFormatter())
    clean_axis(ax)
    panel_label(ax, "(b)")

    ax = axes[2]
    ax.semilogy(
        cycle[callback],
        residual[callback],
        color=ORANGE,
        marker="^",
        markevery=16,
        mfc="white",
        mec=ORANGE,
        mew=0.8,
    )
    ax.axhline(1e-8, color=BLACK, lw=0.9, ls=(0, (3, 2)))
    ax.plot(
        160,
        residual[-1],
        "o",
        color=ORANGE,
        ms=4.6,
        zorder=4,
    )
    ax.set_xlabel("Davidson update")
    ax.set_ylabel(r"$\Vert Hc_m-E_mc_m\Vert\;(E_{\rm h})$")
    ax.set_xlim(0, 160)
    ax.set_ylim(3e-9, 1e-2)
    ax.text(
        4.0,
        1.20e-8,
        r"target $10^{-8}\,E_{\rm h}$",
        ha="left",
        va="bottom",
        color=BLACK,
        fontsize=8.0,
    )
    ax.annotate(
        r"final $6.90\times10^{-7}\,E_{\rm h}$",
        xy=(160.0, residual[-1]),
        xytext=(154.0, 2.5e-4),
        ha="right",
        va="bottom",
        fontsize=8.0,
        color=ORANGE,
        arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": 0.7},
    )
    clean_axis(ax)
    panel_label(ax, "(c)")
    save(fig, "figure1")


def figure_two():
    benchmark_rows = read_csv(DATA / "topk_tracker_basis.csv")
    natural_rows = read_csv(DATA / "topk_natural_basis.csv")
    decomposition_rows = read_csv(
        DATA / "hamiltonian_decomposition_tracker_basis.csv"
    )
    tracker_summary = json.loads(
        (DATA / "compressibility_tracker_basis.json").read_text(
            encoding="utf-8"
        )
    )
    natural_summary = json.loads(
        (DATA / "compressibility_natural_basis.json").read_text(
            encoding="utf-8"
        )
    )
    for label, summary in (
        ("benchmark", tracker_summary),
        ("natural", natural_summary),
    ):
        if not all(bool(value) for value in summary["gates"].values()):
            raise SystemExit(f"{label} summary contains a failed gate")
    if natural_summary["vector_sha256"] != (
        "3855dcb7ce8c44b0f6d9c99fc5787aa2db339002b4d50ec2af8a814c45372d77"
    ):
        raise SystemExit("natural-basis vector hash mismatch")
    gate_topk_csv(benchmark_rows, tracker_summary, "benchmark")
    gate_topk_csv(natural_rows, natural_summary, "natural")
    same_support = json.loads(
        (DATA / "top512_rediagonalization.json").read_text(encoding="utf-8")
    )
    guided_support = json.loads(
        (DATA / "optimized_512_validation.json").read_text(encoding="utf-8")
    )
    natural_same_support = json.loads(
        (DATA / "natural_top512_rediagonalization.json").read_text(
            encoding="utf-8"
        )
    )
    for label, result in (
        ("benchmark same-support", same_support),
        ("natural same-support", natural_same_support),
        ("Hamiltonian-guided", guided_support),
    ):
        if not all(bool(value) for value in result["gates"].values()):
            raise SystemExit(f"{label} 512-determinant control failed a gate")

    def topk_arrays(rows):
        return (
            np.array([int(row["k"]) for row in rows]),
            np.array([float(row["cumulative_norm2"]) for row in rows]),
            np.array(
                [float(row["energy_error_vs_full_eh"]) for row in rows]
            ),
            np.array([float(row["residual_norm_eh"]) for row in rows]),
        )

    k_b, weight_b, energy_b, residual_b = topk_arrays(benchmark_rows)
    k_n, weight_n, energy_n, residual_n = topk_arrays(natural_rows)
    if not np.array_equal(k_b, k_n):
        raise SystemExit("top-K grids differ between orbital representations")
    full_energy = float(tracker_summary["full_state_energy_eh"])
    same_support_error = (
        float(same_support["same_support_optimized_energy_pyscf_eh"])
        - full_energy
    )
    guided_support_error = (
        float(guided_support["energy_pyscf_direct_contraction_eh"])
        - full_energy
    )
    natural_same_support_error = (
        float(natural_same_support["same_support_optimized_energy_pyscf_eh"])
        - full_energy
    )

    k_d = np.array([int(row["k"]) for row in decomposition_rows])
    weight_d = np.array(
        [float(row["retained_norm2"]) for row in decomposition_rows]
    )
    cross_d = np.array(
        [float(row["retained_tail_cross_eh"]) for row in decomposition_rows]
    )
    tail_d = np.array(
        [float(row["tail_h_tail_eh"]) for row in decomposition_rows]
    )
    full_d = np.array(
        [float(row["full_energy_eh"]) for row in decomposition_rows]
    )
    tail_shifted_d = tail_d - (1.0 - weight_d) * full_d
    coupling_contribution = -cross_d / weight_d
    tail_contribution = -tail_shifted_d / weight_d
    decomposed_error = coupling_contribution + tail_contribution
    if np.max(np.abs(decomposed_error - energy_b)) > 6.1e-14:
        raise SystemExit("gauge-invariant energy decomposition failed")

    fig, axes = plt.subplots(2, 2, figsize=(7.05, 4.28))
    fig.subplots_adjust(
        left=0.082, right=0.988, bottom=0.125, top=0.965, hspace=0.43, wspace=0.30
    )

    ax = axes[0, 0]
    ax.semilogx(
        k_b,
        weight_b,
        color=BLUE,
        marker="o",
        mfc="white",
        mec=BLUE,
        mew=0.8,
        label="benchmark orbitals",
    )
    ax.semilogx(
        k_n,
        weight_n,
        color=ORANGE,
        ls=(0, (4, 2)),
        marker="s",
        mfc="white",
        mec=ORANGE,
        mew=0.8,
        label="natural orbitals",
    )
    ax.set_ylabel(r"cumulative norm $W_K$")
    ax.set_xlim(48, 6e6)
    ax.set_ylim(0.0, 1.015)
    ax.legend(loc="lower right", frameon=False, handlelength=2.2)
    clean_axis(ax)
    panel_label(ax, "(a)", x=-0.14)

    ax = axes[0, 1]
    ax.loglog(
        k_b,
        energy_b,
        color=BLUE,
        marker="o",
        mfc="white",
        mec=BLUE,
        mew=0.8,
        label="benchmark orbitals",
    )
    ax.loglog(
        k_n,
        energy_n,
        color=ORANGE,
        ls=(0, (4, 2)),
        marker="s",
        mfc="white",
        mec=ORANGE,
        mew=0.8,
        label="natural orbitals",
    )
    ax.plot(
        512,
        same_support_error,
        marker="D",
        ms=4.7,
        color=BLACK,
        ls="none",
        label="benchmark top-512, relaxed",
        zorder=5,
    )
    ax.plot(
        512,
        guided_support_error,
        marker="*",
        ms=6.4,
        color=GREEN,
        ls="none",
        label="Hamiltonian-guided 512",
        zorder=5,
    )
    ax.plot(
        512,
        natural_same_support_error,
        marker="p",
        ms=5.5,
        color=ORANGE,
        ls="none",
        label="natural top-512, relaxed",
        zorder=5,
    )
    ax.set_ylabel(r"$E_K-E_{\rm Ritz}\;(E_{\rm h})$")
    ax.set_xlim(48, 6e6)
    ax.legend(
        loc="lower left",
        frameon=False,
        handlelength=2.0,
        borderaxespad=0.15,
        fontsize=6.8,
    )
    clean_axis(ax)
    panel_label(ax, "(b)", x=-0.14)

    ax = axes[1, 0]
    ax.loglog(
        k_b,
        residual_b,
        color=BLUE,
        marker="o",
        mfc="white",
        mec=BLUE,
        mew=0.8,
        label="benchmark orbitals",
    )
    ax.loglog(
        k_n,
        residual_n,
        color=ORANGE,
        ls=(0, (4, 2)),
        marker="s",
        mfc="white",
        mec=ORANGE,
        mew=0.8,
        label="natural orbitals",
    )
    ax.set_xlabel(r"retained determinants $K$")
    ax.set_ylabel(r"$\Vert H\phi_K-E_K\phi_K\Vert\;(E_{\rm h})$")
    ax.set_xlim(48, 6e6)
    clean_axis(ax)
    panel_label(ax, "(c)", x=-0.14)

    ax = axes[1, 1]
    ax.semilogx(
        k_d,
        decomposed_error,
        color=BLACK,
        marker="o",
        mfc="white",
        mec=BLACK,
        mew=0.8,
        label=r"total $E_K-E_{\rm Ritz}$",
    )
    ax.semilogx(
        k_d,
        coupling_contribution,
        color=PURPLE,
        ls=(0, (4, 2)),
        marker="^",
        mfc="white",
        mec=PURPLE,
        mew=0.8,
        label="retained--tail coupling",
    )
    ax.semilogx(
        k_d,
        tail_contribution,
        color=GREEN,
        ls=(0, (2, 2)),
        marker="v",
        mfc="white",
        mec=GREEN,
        mew=0.8,
        label="tail quadratic term",
    )
    ax.axhline(0.0, color=GREY, lw=0.65)
    ax.set_xlabel(r"retained determinants $K$")
    ax.set_ylabel(r"contribution to $E_K-E_{\rm Ritz}\;(E_{\rm h})$")
    ax.set_xlim(48, 6e6)
    ax.legend(
        loc="upper right",
        frameon=False,
        handlelength=2.0,
        fontsize=7.2,
    )
    clean_axis(ax)
    panel_label(ax, "(d)", x=-0.14)
    save(fig, "figure2")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    figure_one()
    figure_two()
