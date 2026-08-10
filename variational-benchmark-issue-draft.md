# Quantum Advantage Tracker issue draft — review copy

**Do not submit this file yet.** It mirrors the Tracker's current
“Variational problems submission (existing problem)” form. The proposed issue
reports a new explicit classical variational upper bound for the existing
`2fe_2s_30e_20o` Hamiltonian. It does not claim a new literature minimum,
Davidson convergence, an exact ground-state energy, or quantum advantage.

Proposed issue title:

> [2Fe-2S] (30e, 20o) Davidson Ritz

### Name

[2Fe-2S] (30e, 20o) Davidson Ritz

### Hamiltonian

2fe_2s_30e_20o

### Qubits

_No response_

### Gates

_No response_

### Energy (Eh)

-116.60560912042631

### Low error bound (Eh)

_No response_

### High error bound (Eh)

_No response_

### Method

Checkpointed full-space Davidson Ritz (PySCF 2.14.0)

### Method proof

Reproduction repository:
https://github.com/Osgood001/2fe2s-variational-benchmark

This result is a Rayleigh quotient of an explicit CI vector in the complete
CAS(30e,20o), fixed-`M_S=0` determinant space with
`(N_alpha,N_beta)=(15,15)`. The Hamiltonian dimension is `240,374,016`.
PySCF's matrix-free `direct_spin1` Davidson solver applies the full
Slater–Condon Hamiltonian without forming the matrix.

The calculation was regenerated from the public FCIDUMP, not from a contestant
or literature wavefunction. The deterministic initial CI tensor had two
nonzero entries, `c[0,0]=1.00001` and
`c[15503,15503]=-1e-5`. A 110-update Davidson run generated the first saved
vector. Its hash-identified endpoint was then continued for 2, 16, and 160
updates, respectively, with a wider 64-vector Davidson restart space. The
stage endpoint energies were:

- A: `-116.6055095343891 Eh`;
- B1: `-116.6055107620259 Eh`;
- B2: `-116.6055790186671 Eh`;
- C: `-116.60560912042631 Eh`.

The final numerical evidence is:

- FCIDUMP SHA-256:
  `bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7`;
- final vector SHA-256:
  `45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945`;
- vector shape/dtype: `15504 x 15504`, little-endian float64
  (`240,374,016` real coefficients; `1,922,992,256` bytes);
- squared norm: `0.999999999999999`;
- independently contracted Rayleigh energy:
  `-116.60560912042631 Eh`;
- residual norm `||Hc-Ec||`: `6.901557706123275e-7 Eh`;
- energy variance: `4.763149877094956e-13 Eh^2`;
- independent RDM energy: `-116.60560912042571 Eh`;
- absolute RDM/direct energy difference: `5.97e-13 Eh`;
- `<S^2>`: `2.8569746519844883e-9`;
- PySCF solver convergence flag: `false`.

The energy was recomputed by a fresh Hamiltonian action on the saved normalized
vector; it was not copied from a solver log. The nonzero residual is reported
as a convergence diagnostic and is not converted into an unsupported
two-sided error bar. The direct Rayleigh quotient remains a variational upper
bound in the stated particle-number and spin-projection sector.

The exact final vector is identified above by its SHA-256. It is not committed
to the reproduction repository because it is 1.79 GiB; a separate durable
archive can be added during verification if requested.

### Quantum runtime (seconds)

_No response_

### Classical runtime (seconds)

32746

### Compute resources (quantum)

_No response_

### Compute resources (classical)

64 CPU threads; 512 GiB memory allocation; PySCF 2.14.0

### Notes

Explicit unconverged Ritz vector; fresh-action replay

### Authors

Shigang Ou

### Institutions

Institute of Physics, Chinese Academy of Sciences

---

## Reviewer context — do not paste into the Tracker form

The current verified Tracker entry is HCI energy `-116.605425 Eh` in
[issue #187](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/187).
The proposed value is lower by `0.00018412042631 Eh` (approximately
`184.12 microEh`). This is why the update is useful to the Tracker.

It is not a new literature-low result. It agrees, at the reported precision,
with the DMRG value `-116.6056091 Eh` published by Li and Chan in 2017
([DOI: 10.1021/acs.jctc.7b00270](https://doi.org/10.1021/acs.jctc.7b00270)).
The narrow claim is therefore “lower than the value currently represented in
the Tracker,” not “lower than every published classical calculation.”

The `32746 s` runtime is the rounded sum of the recorded end-to-end wall times
for A (`7667.02 s`), B1 (`258.25 s`), B2 (`1326.97 s`), and C
(`23494.16 s`). It describes the exact staged construction, not an optimized
time-to-solution benchmark.

## Pre-submission checklist — do not paste into the Tracker form

- [ ] Review and approve the scientific wording.
- [ ] Deposit the exact 1.79-GiB final vector in a durable, unauthenticated
      archive if the Tracker maintainers require it for verification.
- [ ] Confirm that the downloaded archive vector matches SHA-256
      `45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945`.
- [ ] Recheck the live Tracker incumbent immediately before submission.
- [ ] Keep low/high error fields empty unless a rigorous two-sided bound is
      established.
- [ ] Do not describe the state as exact, converged, or a literature SOTA.
