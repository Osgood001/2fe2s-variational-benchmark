# Independent referee report — round 1

## Decision

**REJECT / NOT PUBLISHABLE.** The numerical object appears credible, but the
submitted PDF is an internal computational audit rather than a Physical Review
Letters manuscript. A scientific and graphical rewrite is required.

## Evidence accepted by the referee

- The official and working FCIDUMP bytes are identical.
- The explicit Ritz vector is normalized; direct Rayleigh contraction and RDM
  reconstruction differ by approximately \(6\times10^{-13}E_h\).
- The residual is \(6.90\times10^{-7}E_h\), and the spin diagnostics are
  consistent with the target singlet.
- The orbital-rotation equivalence of the two public Hamiltonian
  representations is well supported.
- The improvement of approximately \(184.1\,\mu E_h\) over the current Tracker
  HCI entry is numerically credible.

## Major comments

1. Establish a broad scientific question: why the CAS(30e,20o) [2Fe--2S]
   Hamiltonian is a quantum-simulation stress test, what is missing from
   existing classical baselines, and what an explicit determinant-space Ritz
   state changes.
2. Agreement with the rounded literature DMRG energy is not by itself a new
   PRL-level physical result. Add a generalizable wave-function compressibility
   analysis: top-\(K\) accumulated norm, directly recomputed top-\(K\) energy
   and residual, and natural-orbital occupations. Determine whether coefficient
   weight concentrates much faster than correlation energy.
3. Consistently say "the full configuration-interaction space within
   CAS(30e,20o), in the fixed \(M_S=0\) sector" and "Ritz variational upper
   bound." Never imply full-molecule/full-basis FCI or exact convergence.
4. Report the HCI improvement only to the incumbent's precision. Treat the
   literature DMRG value as a rounding interval; do not claim a meaningful
   \(0.020\,\mu E_h\) improvement or a two-sided error bar.
5. Call the second PySCF calculation an "independent replay" or "independent
   execution," not an independent software implementation.
6. Describe the RDM1 tolerance change as a post-hoc numerical-stability study,
   not a predeclared gate. Preserve the initial failure and its diagnostics in
   Supplemental Material if it remains relevant.
7. State the orbital transformation equations and numerical tensor residuals;
   delete the low-information two-point threshold plot.
8. Remove operational material from the Letter: attempt/Job identifiers, CLI
   and platform failures, dependency history, hash tables, submission advice,
   and upload instructions.
9. Resolve authors, affiliations, contribution approval, and an unauthenticated
   permanent data archive before actual journal submission.

## Figure requirements

- Delete the workflow diagram and two-point threshold graphic.
- Correct the duplicated \(10^6\) axis conversion.
- Replace energy bars by point/interval graphics; represent the DMRG reference
  as a rounding band.
- Use at most two main composite figures: one for the result/convergence and
  one for wave-function compressibility/strong correlation.
- Use English panel labels, at least 8 pt final lettering, color-blind-safe
  colors plus redundant markers/line styles, no large background grids, and
  vector outputs generated from frozen data and scripts.

## Publishability gates

- [ ] Title/abstract precisely limit the claim to a CAS(30e,20o) Ritz upper bound.
- [ ] No quantum-advantage, exact-FCI, 512-determinant SOTA, or DMRG-beating claim.
- [ ] A new, generalizable compressibility/strong-correlation result is present.
- [ ] Hamiltonian source, energy convention, and orbital-equivalence equations are explicit.
- [ ] Precision, DMRG rounding, nonconvergence, residual, replay independence, and post-hoc RDM tolerance are represented honestly.
- [ ] Main text is English, physics-led, four-page PRL length, and free of internal operations.
- [ ] Main figures satisfy the point/interval, convergence, compressibility, occupation, readability, and reproducibility requirements.
- [ ] Public data/code archive, author list, affiliations, contribution approval, and current incumbent check are complete before external submission.

The referee will only return **PUBLISHABLE** after every scientific,
presentation, and reproducibility gate is satisfied.
