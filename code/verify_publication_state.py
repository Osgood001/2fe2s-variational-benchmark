#!/usr/bin/env python3
"""Fresh-process verification for the explicit 2Fe-2S Ritz state.

This script deliberately does not import the continuation driver.  It reloads
the FCIDUMP and saved vector, rebuilds the Hamiltonian contraction, and records
the direct, reduced-density-matrix, residual, norm, and spin observables used
in the manuscript.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np
import pyscf
from pyscf import ao2mo, fci, lib
from pyscf.fci import spin_op
from pyscf.tools import fcidump


EXPECTED_FCIDUMP_SHA256 = (
    "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7"
)
PUBLIC_HCI_BOUND_EH = -116.605425


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


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: verify_publication_state.py FCIDUMP VECTOR_NPY OUTPUT_JSON"
        )

    started = time.time()
    fcidump_path = Path(sys.argv[1]).resolve()
    vector_path = Path(sys.argv[2]).resolve()
    output_path = Path(sys.argv[3]).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    threads = int(os.environ.get("OMP_NUM_THREADS", "64"))
    lib.num_threads(threads)

    fcidump_sha256 = sha256_file(fcidump_path)
    vector_sha256_before = sha256_file(vector_path)
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

    shape = (math.comb(norb, nelec[0]), math.comb(norb, nelec[1]))
    civec = np.load(vector_path, allow_pickle=False)
    if civec.shape != shape or civec.dtype != np.float64:
        raise RuntimeError(
            f"vector shape/dtype mismatch: {civec.shape} {civec.dtype}"
        )

    norm2 = float(np.vdot(civec.ravel(), civec.ravel()).real)
    hc = np.asarray(
        fci.direct_spin1.contract_2e(h2e, civec, norb, nelec), dtype=np.float64
    ).reshape(shape)
    electronic_energy = float(np.vdot(civec.ravel(), hc.ravel()).real / norm2)
    energy_eh = electronic_energy + ecore
    hc -= electronic_energy * civec
    residual_norm_eh = float(np.linalg.norm(hc.ravel()) / math.sqrt(norm2))
    del hc

    # For N_alpha=N_beta, a pure even-spin state has C_ab = C_ba.  Evaluate
    # this without materializing another full 1.79-GiB array.
    exchange_overlap = float(np.einsum("ij,ji->", civec, civec, optimize=True))
    transpose_defect2 = max(0.0, 2.0 * (norm2 - exchange_overlap) / norm2)
    transpose_relative_defect = math.sqrt(transpose_defect2)

    # S^2 is evaluated by an independent operator contraction, not inferred
    # solely from exchange symmetry.
    spin_square, _multiplicity_raw = spin_op.spin_square0(civec, norb, nelec)
    spin_square = float(spin_square / norm2)
    multiplicity = float(2.0 * math.sqrt(spin_square + 0.25))

    # Reconstruct the energy through the independent reduced-density-matrix
    # kernels.  This is deliberately a different code path from contract_2e
    # and catches integral-index/factor mistakes that a repeated contraction
    # could share with the continuation driver.
    dm1_rdm, dm2_rdm = fci.direct_spin1.make_rdm12(civec, norb, nelec)
    dm1_rdm /= norm2
    dm2_rdm /= norm2
    eri_full = ao2mo.restore(1, eri, norb)
    rdm_electronic_energy = float(
        np.einsum("pq,qp->", h1, dm1_rdm, optimize=True)
        + 0.5 * np.einsum("pqrs,pqrs->", eri_full, dm2_rdm, optimize=True)
    )
    rdm_energy_eh = rdm_electronic_energy + ecore
    rdm_vs_contraction_abs_eh = abs(rdm_energy_eh - energy_eh)
    del dm2_rdm, eri_full

    dm1a, dm1b = fci.direct_spin1.make_rdm1s(civec, norb, nelec)
    dm1a /= norm2
    dm1b /= norm2
    dm1 = dm1a + dm1b

    # make_rdm12 and make_rdm1s use distinct compiled accumulation paths.  On
    # the 240,374,016-coefficient state their elementwise difference is small
    # but no longer at the four-orbital preflight's machine-epsilon scale.
    # Preserve enough diagnostics to distinguish accumulation-order roundoff
    # from a transpose, symmetry, trace, or normalization error.
    dm1_difference = dm1 - dm1_rdm
    dm1_difference_abs = np.abs(dm1_difference)
    rdm1_path_max_index = tuple(
        int(i) for i in np.unravel_index(np.argmax(dm1_difference_abs), dm1.shape)
    )
    rdm1_path_max_abs_difference = float(dm1_difference_abs[rdm1_path_max_index])
    rdm1_path_signed_difference_at_max = float(
        dm1_difference[rdm1_path_max_index]
    )
    rdm1_path_value_spin_resolved_at_max = float(dm1[rdm1_path_max_index])
    rdm1_path_value_rdm12_at_max = float(dm1_rdm[rdm1_path_max_index])
    rdm1_path_frobenius_difference = float(np.linalg.norm(dm1_difference))
    rdm1_path_relative_frobenius_difference = float(
        rdm1_path_frobenius_difference
        / max(float(np.linalg.norm(dm1)), float(np.linalg.norm(dm1_rdm)))
    )
    rdm1_path_relative_max_difference = float(
        rdm1_path_max_abs_difference
        / max(float(np.max(np.abs(dm1))), float(np.max(np.abs(dm1_rdm))))
    )
    rdm1_spin_resolved_symmetry_max_abs = float(np.max(np.abs(dm1 - dm1.T)))
    rdm1_rdm12_symmetry_max_abs = float(np.max(np.abs(dm1_rdm - dm1_rdm.T)))
    rdm1_spin_resolved_trace = float(np.trace(dm1))
    rdm1_rdm12_trace = float(np.trace(dm1_rdm))
    rdm1_trace_path_abs_difference = abs(
        rdm1_spin_resolved_trace - rdm1_rdm12_trace
    )
    rdm1_one_body_energy_spin_resolved_eh = float(
        np.einsum("pq,qp->", h1, dm1, optimize=True)
    )
    rdm1_one_body_energy_rdm12_eh = float(
        np.einsum("pq,qp->", h1, dm1_rdm, optimize=True)
    )
    rdm1_one_body_energy_path_abs_difference_eh = abs(
        rdm1_one_body_energy_spin_resolved_eh - rdm1_one_body_energy_rdm12_eh
    )
    natural_occupations = np.linalg.eigvalsh((dm1 + dm1.T) * 0.5)[::-1]
    spin_density_frobenius = float(np.linalg.norm(dm1a - dm1b))

    # Blocked coefficient statistics avoid a full-sized temporary array.
    sum_weight2 = 0.0
    shannon_nats = 0.0
    max_abs_coefficient = 0.0
    counts = {"1e-2": 0, "1e-3": 0, "1e-4": 0, "1e-5": 0, "1e-6": 0}
    for row0 in range(0, shape[0], 128):
        block = civec[row0 : row0 + 128]
        weights = np.square(block, dtype=np.float64) / norm2
        sum_weight2 += float(np.sum(weights * weights))
        nonzero = weights[weights > 0.0]
        shannon_nats -= float(np.sum(nonzero * np.log(nonzero)))
        max_abs_coefficient = max(max_abs_coefficient, float(np.max(np.abs(block))))
        abs_block = np.abs(block)
        for label, threshold in (
            ("1e-2", 1e-2),
            ("1e-3", 1e-3),
            ("1e-4", 1e-4),
            ("1e-5", 1e-5),
            ("1e-6", 1e-6),
        ):
            counts[label] += int(np.count_nonzero(abs_block >= threshold))

    vector_sha256_after = sha256_file(vector_path)
    if vector_sha256_after != vector_sha256_before:
        raise RuntimeError("vector file changed during verification")

    result = {
        "schema_version": "2fe2s-publication-state-verification/v2",
        "generated_at_unix": time.time(),
        "host": platform.node(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "pyscf_version": pyscf.__version__,
        "threads": lib.num_threads(),
        "fcidump_path": str(fcidump_path),
        "fcidump_sha256": fcidump_sha256,
        "vector_path": str(vector_path),
        "vector_bytes": vector_path.stat().st_size,
        "vector_sha256": vector_sha256_after,
        "norb": norb,
        "nelec": list(nelec),
        "ms2": ms2,
        "dimension": shape[0] * shape[1],
        "vector_shape": list(shape),
        "norm2": norm2,
        "ecore_eh": ecore,
        "electronic_energy_eh": electronic_energy,
        "energy_rayleigh_eh": energy_eh,
        "residual_norm_eh": residual_norm_eh,
        "energy_variance_eh2": residual_norm_eh**2,
        "rdm_electronic_energy_eh": rdm_electronic_energy,
        "rdm_energy_eh": rdm_energy_eh,
        "rdm_vs_contraction_abs_eh": rdm_vs_contraction_abs_eh,
        "rdm1_path_max_abs_difference": rdm1_path_max_abs_difference,
        "rdm1_path_max_index": list(rdm1_path_max_index),
        "rdm1_path_signed_difference_at_max": rdm1_path_signed_difference_at_max,
        "rdm1_path_value_spin_resolved_at_max": rdm1_path_value_spin_resolved_at_max,
        "rdm1_path_value_rdm12_at_max": rdm1_path_value_rdm12_at_max,
        "rdm1_path_frobenius_difference": rdm1_path_frobenius_difference,
        "rdm1_path_relative_frobenius_difference": (
            rdm1_path_relative_frobenius_difference
        ),
        "rdm1_path_relative_max_difference": rdm1_path_relative_max_difference,
        "rdm1_spin_resolved_symmetry_max_abs": (
            rdm1_spin_resolved_symmetry_max_abs
        ),
        "rdm1_rdm12_symmetry_max_abs": rdm1_rdm12_symmetry_max_abs,
        "rdm1_spin_resolved_trace": rdm1_spin_resolved_trace,
        "rdm1_rdm12_trace": rdm1_rdm12_trace,
        "rdm1_trace_path_abs_difference": rdm1_trace_path_abs_difference,
        "rdm1_one_body_energy_spin_resolved_eh": (
            rdm1_one_body_energy_spin_resolved_eh
        ),
        "rdm1_one_body_energy_rdm12_eh": rdm1_one_body_energy_rdm12_eh,
        "rdm1_one_body_energy_path_abs_difference_eh": (
            rdm1_one_body_energy_path_abs_difference_eh
        ),
        "improvement_vs_public_hci_eh": PUBLIC_HCI_BOUND_EH - energy_eh,
        "exchange_overlap": exchange_overlap / norm2,
        "transpose_relative_defect": transpose_relative_defect,
        "spin_square": spin_square,
        "spin_multiplicity": multiplicity,
        "rdm1_trace": rdm1_spin_resolved_trace,
        "spin_density_frobenius": spin_density_frobenius,
        "natural_occupations_descending": natural_occupations.tolist(),
        "max_abs_coefficient": max_abs_coefficient,
        "inverse_participation_ratio": sum_weight2,
        "effective_determinant_count_ipr": 1.0 / sum_weight2,
        "coefficient_shannon_entropy_nats": shannon_nats,
        "coefficient_effective_count_shannon": math.exp(shannon_nats),
        "coefficient_counts_abs_ge": counts,
        "wall_s": time.time() - started,
        "process_max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    atomic_json(output_path, result)
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
