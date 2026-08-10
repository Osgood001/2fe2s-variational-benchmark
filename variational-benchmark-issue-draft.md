# Variational-problem issue draft

This file mirrors the fields of the existing-problem submission form.  Replace
the three marked attribution/archive fields before submitting.

## Name

[2Fe-2S] (30e, 20o) explicit Davidson Ritz variational bound

## Hamiltonian

2fe_2s_30e_20o

## Qubits

No response

## Gates

No response

## Energy (Eh)

-116.60560912042631

## Low error bound (Eh)

No response

## High error bound (Eh)

No response

## Method

Checkpointed Davidson Ritz minimization in the complete CAS(30e,20o),
fixed-M_S=0 determinant space (PySCF 2.14.0)

## Method proof

We minimize the Rayleigh quotient in the complete determinant space of the
specified CAS(30e,20o) Hamiltonian, restricted to
`(N_alpha,N_beta)=(15,15)`.  The matrix-free `direct_spin1` Davidson procedure
produces explicit Ritz vectors; the reported energy was recomputed by applying
the Hamiltonian to the final saved vector rather than copied from a solver log.

The frozen numerical evidence is:

- FCIDUMP SHA-256:
  `bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7`;
- vector SHA-256:
  `45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945`;
- vector shape and dtype: `15504 x 15504`, little-endian float64
  (`240,374,016` real coefficients);
- squared norm: `0.999999999999999`;
- directly contracted Rayleigh energy: `-116.60560912042631 Eh`;
- residual norm `||Hc-Ec||`: `6.901557706123275e-7 Eh`;
- energy variance: `4.763149877094956e-13 Eh^2`;
- reduced-density-matrix energy: `-116.60560912042571 Eh`;
- absolute RDM/direct energy difference: `5.9686e-13 Eh`;
- `<S^2>`: `2.8569746519844883e-9`;
- solver convergence flag: `false`.

The direct Rayleigh quotient of this normalized state is a variational upper
bound for the stated Hamiltonian and symmetry sector.  The nonzero residual is
reported as a convergence diagnostic; it is not converted into an unsupported
two-sided energy error.  The Supplemental Material documents the complete
staged calculation, Hamiltonian-basis equivalence, fresh-process replay, and
reproduction commands.

Manuscript, scripts, compact data, and final-vector archive:
`{{PUBLIC_ARTIFACT_URL}}`

## Quantum runtime (seconds)

No response

## Classical runtime (seconds)

No response

## Compute resources (quantum)

No response

## Compute resources (classical)

No response

## Notes

Explicit CAS Ritz vector; residual replayed in a fresh process.  On 2026-08-09
the current verified Tracker entry remained the HCI value
`-116.605425 Eh` in
<https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/187>.
The new Rayleigh value is therefore approximately `184 microEh` lower at the
precision supported by that entry.  It agrees with, and is not claimed to
improve upon, the published rounded DMRG value `-116.6056091 Eh`.

## Authors

{{APPROVED_AUTHORS}}

## Institutions

{{APPROVED_INSTITUTIONS}}

---

## Pre-submission checks (remove this section from the submitted issue)

- [ ] Replace all three `{{...}}` fields.
- [ ] Obtain contribution and attribution approval from every contributor.
- [ ] Deposit the final vector, FCIDUMP, scripts, environment locks, compact
      outputs, Letter, and Supplemental Material in a durable archive.
- [ ] Verify the archive URL without authentication and record its checksum.
- [ ] Recheck the public incumbent immediately before submission.
- [ ] Keep low/high error fields empty unless a rigorous two-sided bound is
      established.
- [ ] Do not describe the state as exact or converged.
