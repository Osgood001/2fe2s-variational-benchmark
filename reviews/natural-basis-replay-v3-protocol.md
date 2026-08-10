# Independent protocol review — natural-basis replay v3

## Initial decision: NO-GO

The referee identified two blockers before the replay could be submitted:

1. The frozen primary-basis occupation CSV lacked a preflight SHA-256 gate.
2. RDM symmetry and agreement of occupation eigenvalues did not establish that
   the transformed CI tensor was expressed in the stated natural-orbital basis.

The referee accepted revising the occupation-trace absolute tolerance from
`1e-8` to `2e-8` only if the change was disclosed as a response to the
preliminary-run diagnostic, frozen before v3, and never adjusted after seeing
v3.

## Revision

- The wrapper and analyzer now check the occupation CSV before any FCIDUMP or
  full-vector analysis. Its frozen SHA-256 is
  `2e9632ff525337eae6329eac5cd5e39b71bd08af85291a335033296b5a732b7d`.
- The summary records the expected and observed hashes and a hash gate.
- The replay gates the raw RDM symmetry defect, the maximum off-diagonal
  element of the symmetrized RDM, its complete elementwise difference from
  `diag(frozen_primary_occupations)`, and the occupation-eigenvalue difference.
- Both the trace and RDM elementwise tolerances were frozen at `2e-8` before
  submission.
- Python 3.12 compilation, shell syntax, and the remote reference-file hash
  were checked before submission.

## Resubmission decision: GO

The same referee returned an exact `GO` after inspecting the revised scripts.
The scientific thresholds are frozen for Bohrium job `20529398`.
