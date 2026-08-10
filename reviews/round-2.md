# Independent referee report — round 2

## Decision

**NOT PUBLISHABLE.**

The revision is dramatically better and is now a conventional scientific
manuscript rather than an audit report.  No blocking numerical error was found
in the Ritz energy, replay, Hamiltonian equivalence, or top-\(K\) contractions,
and the graphics meet the publication-quality threshold.  Two scientific
issues still block a PRL-level recommendation: the determinant-compressibility
result is basis dependent but is presented too broadly, and its algorithmic
implication lacks a Hamiltonian-aware control.  The complete vector also lacks
a permanent public archive.

## Round-1 disposition

| Round-1 requirement | Disposition | Round-2 assessment |
|---|---|---|
| Physics-led question rather than competition history | Satisfied | The introduction now motivates strong correlation, classical baselines, and representation complexity. |
| New compressibility/strong-correlation result | Partially satisfied | Exact top-\(K\) contractions and occupations are substantial new evidence, but only in one orbital basis and without an optimized-\(K\) control. |
| Limit FCI claim to CAS(30e,20o), \(M_S=0\) | Satisfied | Scope is explicit in Letter and Supplemental Material. |
| Disclose nonconvergence and residual | Satisfied | `converged=false`, the \(10^{-8}E_h\) target, and the \(6.90\times10^{-7}E_h\) residual are all prominent. |
| Do not claim a certified two-sided error | Satisfied | The text correctly identifies a Rayleigh upper bound only. |
| Treat DMRG as a rounding interval | Satisfied in prose; minor figure issue | The prose is correct, but plotted `+/-0.05` still visually resembles an uncertainty estimate. |
| Restrict HCI comparison to available precision | Partially satisfied | Prose acknowledges the limitation, but `184.1` and `-184.12` overstate precision relative to the six-decimal HCI entry. |
| Say fresh-process replay, not independent implementation | Satisfied | The supplement explicitly says both paths use PySCF. |
| Present RDM tolerance change as post-hoc stability study | Satisfied | The initial failure and expanded diagnostics are transparently retained. |
| Give orbital-equivalence equations and tensor residuals | Satisfied | Equations and complete-tensor RMS/max errors are clear and reproducible. |
| Remove operational narrative from Letter | Satisfied | Job IDs, dependencies, hashes, and restart details are confined to the supplement/bundle. |
| Replace workflow/threshold graphics | Satisfied | Both were removed. |
| Correct energy-axis units | Satisfied | No duplicated conversion remains. |
| Use point/interval energy comparison | Satisfied with minor labeling revision | Fig. 1(a) is conceptually correct. |
| Add optimized \(K=512\) comparison | Not satisfied | The value exists in frozen data but is absent from the paper and figure. |
| Two scientific composite figures | Satisfied | Figures are data-driven and materially support the argument. |
| At least 8-pt final lettering and redundant encoding | Satisfied | Current figures are readable, vector-based, color-blind compatible, and use markers/line styles. |
| Frozen data and figure-generation scripts | Satisfied | The figure pipeline maps cleanly to CSV/JSON. |
| Permanent public archive, including full vector | Not satisfied | Remote storage and “available from the authors” are not a permanent, unauthenticated archive. |
| Authors and affiliations resolved | Administrative pending | Acceptable for anonymous scientific review, but required before submission. |
| Current incumbent checked | Satisfied | The official data still list HCI at \(-116.605425E_h\). |

## Numerical correctness accepted

- The streaming procedure returns the exact global largest squared amplitudes.
- Every projected tensor is explicitly populated, normalized, and contracted.
- The overlap equals \(\sqrt{W_K}\).
- At \(K=4{,}194{,}304\), \(W_K=0.9856778276\), overlap is
  0.992813088, energy error is 0.0177050366 \(E_h\), and residual is
  0.1865470 \(E_h\).
- Natural occupations have the required trace.
- Direct and RDM energies agree to \(5.97\times10^{-13}E_h\).
- The corrected \(H_{\rm el}\) notation removes the core-energy ambiguity.
- Orbital-equivalence evidence remains convincing.

## Major comments

### 1. Test orbital-basis dependence

Top-\(K\) determinant compressibility is not invariant under active-active
orbital rotations.  Repeat the analysis in at least one meaningful alternate
representation (the Li--Chan/QMB basis already used in the equivalence test, or
preferably natural orbitals), or supply a comparably strong basis-sensitivity
analysis.  Plot both curves.  If the separation persists, the result becomes
substantially stronger; if it changes, orbital dependence should become the
result.  Merely narrowing every claim to the fixed basis would make a useful
benchmark note, not the current PRL-level claim.

### 2. Add a Hamiltonian-aware control

The current experiment proves that preserving maximum overlap does not preserve
energy, but it does not show what coefficient reoptimization or Hamiltonian
selection repairs.  At minimum, add the independently optimized 512-determinant
state already frozen in the data:

\[
E_{\mathrm{top-}512}=-116.0003434451E_h,\qquad
E_{\mathrm{optimized}\ 512}=-116.3704626758E_h.
\]

The control recovers about 370 m\(E_h\) beyond amplitude projection but remains
about 235 m\(E_h\) above the full Ritz value.  Put it in Fig. 2(b), document
independent provenance, and distinguish fixed-coefficient compression, support
selection plus coefficient reoptimization, and adaptive SCI.  Preferably also
rediagonalize within several top-\(K\) supports.

### 3. Demonstrate or remove the mechanism

For \(p=P_Kc\) and \(q=(1-P_K)c\), directly report
\(\langle p|H|p\rangle\), \(\langle q|H|q\rangle\), and
\(2\operatorname{Re}\langle p|H|q\rangle\) for representative \(K\).  This
will distinguish coherent retained--tail coupling, tail self-energy, and
normalization.  Otherwise remove the off-diagonal-coherence mechanism language.

### 4. PRL significance

The combination of basis-robust separation, a Hamiltonian-aware control, and a
cross-term mechanism analysis could make the benchmark consequence
Letter-worthy.  A single-basis fidelity-versus-energy observation is not yet
sufficient.  The connection to amplitude ranking in the cited HI-VQE work must
state exactly which step is diagnosed and which rediagonalization steps are not
modeled.

### 5. Archive the scientific object

Deposit the 1.79-GiB vector and full reproduction bundle at a permanent,
unauthenticated accession.  Include filename, size, dtype, shape, SHA-256,
FCIDUMP source revision, code, compact data, license, citation instructions,
and a complete manifest.

## Minor comments

1. Round the HCI improvement to approximately 184 microhartree; remove false
   digits in Fig. 1(a).
2. Render DMRG as a labeled published rounding interval, not `+/-` uncertainty.
3. Confirm the HCI number is variational before calling it an upper bound;
   otherwise say entry/value.
4. Replace Fig. 2(b)'s “energy scale” by an explicit ordinate.
5. Define or remove the 0.5--1.5 occupation shading criterion.
6. Say “tail relevant to the energy,” not generic observables.
7. Say “complete staged protocol,” not deterministic schedule.
8. Recheck the rendered 600-character abstract limit.
9. Supply a retrievable source/version for the second FCIDUMP in the archive.
10. Populate PDF title/author metadata in the submission build.

## Exact gates for round 3

- [ ] Multi-basis or comparably strong basis-sensitivity analysis.
- [ ] Independently validated Hamiltonian-aware \(K=512\) control; preferably
      rediagonalized energies in several top-\(K\) supports.
- [ ] Retained/tail/cross Hamiltonian decomposition or removal of mechanism.
- [ ] Precisely scoped SCI, sampling, and compression claims.
- [ ] Honest HCI/DMRG precision.
- [ ] Permanent public accession for the vector and bundle.
- [ ] Minor labels, abstract, source retrieval, and metadata corrected.
- [ ] Authors, affiliations, contribution approval, and authorization resolved
      before external submission.
