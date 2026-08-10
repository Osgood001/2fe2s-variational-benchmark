# Final publication review — OHEC log

## Inquiry

- Question: Are the frozen final Letter, Supplemental Material, figures, data, code, and issue draft scientifically and reproducibly publishable at PRL level?
- Scope and constraints: Read-only referee review of the exact final artifacts and Bohrium job 20529398; anonymous author/accession placeholders are permitted, but any scientific or submission blocker must be reported with exact evidence.
- Success or stopping criterion: Every material claim is consistent with frozen numerical evidence and definitions; variational and literature comparisons are bounded honestly; v2/v3 tolerance history is transparent; all reproduction links resolve; no substantive blocker remains.

## Evidence index

| ID | Source or artifact | What it can establish | Limitations |
|---|---|---|---|
| E1 | `main.tex`, `supplement.tex`, `results.tex`, compiled PDF | Exact scientific claims, definitions, figures, and presentation | Does not independently validate numbers |
| E2 | `bohr-results/20529398/` and `data/` | Frozen v3 numerical outcomes and gates | Numerical evidence remains bounded to the stated active-space Hamiltonian |
| E3 | `analysis/`, `code/`, `environment/`, `README.md` | Reproduction logic, provenance, and executable definitions | External vector availability depends on the pending accession |
| E4 | `variational-benchmark-issue-draft.md`, official tracker/literature references | Scope and accuracy of the proposed public issue | Public submission metadata remain intentionally pending |

## Loop 001

### Observation

- Evidence: The parent reports that v3 job 20529398 exited zero with all 18 checks passing and requests an artifact-level verification rather than reliance on that report.
- Source/provenance: Parent task message; user-reported evidence pending direct inspection.
- Direct observation vs interpretation: This is a reported observation, not yet independently verified.

### Hypothesis

- H1: The revised artifact closes the basis-dependence, Hamiltonian-aware-control, tolerance-transparency, and reproducibility blockers from rounds 1–2.
- Prediction if true: Frozen outputs match all manuscript numbers; the Letter distinguishes projected top-K, same-support rediagonalization, and separately optimized supports; literature values are quoted at justified precision; the v2 failure and v3 tolerance are disclosed; all archive paths resolve except authorized placeholders.
- Falsifier: Any unsupported variational/SOTA claim, mismatched number or caption, hidden post-hoc tolerance, missing frozen input/hash, or reproduction command that cannot generate a claimed result.
- Competing hypotheses: The numerics are correct but the manuscript overstates their generality; or the manuscript is scientifically sound but the reproduction bundle is incomplete.

### Experiment

- Design: Cross-walk every quantitative claim and figure to compact output records and code definitions, then inspect the rendered PDF and all reproduction entry points.
- Control and changed variable: Compare tracker-orbital and natural-orbital representations of the same hash-locked Ritz state; distinguish fixed-support coefficient changes from support changes.
- Decision rule declared in advance: Return `PUBLISHABLE` only if no substantive scientific or reproducibility blocker survives; otherwise enumerate exact blockers.

### Result

- Raw result: Pending systematic inspection.
- Artifact or conversation reference: E1–E4.
- Confounders: Author identities and permanent accession are intentionally withheld pending authorization and will be assessed only as submission logistics.

### Conclusion

- Status: inconclusive
- Confidence and evidence boundary: No verdict before direct artifact inspection.
- Remaining uncertainty: Claim/data consistency, literature precision, graphical consistency, and complete replay closure.
- Next experiment: Parse the manuscript and frozen outputs, then visually inspect the PDF.

### Reflection

The strongest temptation is to infer publishability from the successful v3 gate. The final review must instead test whether the prose makes claims no numerical gate was designed to establish.

## Loop 002

### Observation

- The exact final combined PDF has SHA-256 `8118fa1ab087b1ce4fb07093797f3c2df27ca80dae612f1472d9f2ff1249f368`, is ten pages, embeds its fonts, and renders the two main figures legibly.
- The frozen job 20529398 completed successfully and its compact summary reports all 18 gates as passing. The numerical records in `data/` are hash-identical to the corresponding job outputs. The manuscript's final energy, residual, natural-basis invariance checks, top-​K projections, same-support rediagonalizations, guided support, natural top-512 control, occupation tolerance, and literature comparisons agree with those records at the quoted precision.
- The Letter accurately confines its variational statement to CAS(30e,20o) and the fixed spin sector, reports `converged=false` through the residual discussion rather than claiming an exact eigenvalue, and explicitly says that the rounded DMRG result is matched rather than improved.
- Three final-text problems remain: (i) the rendered abstract is about 647 characters including spaces, above PRL's explicit 600-character limit; (ii) `supplement.tex:242-243` sends the reader to Eq. (2) of the Letter for the orbital-integral transformation, but Letter Eq. (2) is the Rayleigh quotient/residual definition at `main.tex:93-100`; and (iii) the decomposition implementation does not independently contract the tail as the prose claims. In `analysis/analyze_state_compressibility.py:303-309`, `tail_term` is defined by subtracting the retained and cross terms from the full energy, which makes `decomposition_closure` algebraic by construction. Reconstructing the shifted Eq. (5) closure from the frozen tracker CSV gives a maximum `6.039613253960852e-14 Eh`, not the claimed `2.3e-14 Eh` across both bases.
- The revised data-availability wording now accurately distinguishes locally bundled derived natural-orbital artifacts from immutable public locations and hashes of the two source FCIDUMPs. The remaining author/accession placeholders are declared openly and are logistical, not a defect in the frozen scientific analysis.

### Hypothesis

- H2: The numerical study is scientifically sound, but the exact final submission is not yet publishable because one mandatory PRL format rule and two checkable method-description statements are wrong.
- Prediction if true: Correcting these three local text issues will not change any energy, score, figure, or scientific conclusion, after which no substantive blocker should remain.
- Falsifier: A compliant abstract count, a correct Letter Eq. (2) transformation, or an archived independent `Hq` contraction supporting the stated closure.
- Competing hypotheses: The decomposition sentence could be treated as harmless shorthand, but its use of “Direct contractions verify” is stronger than the archived computation and is therefore not acceptable in a reproducibility-focused final manuscript.

### Experiment

- Compared the rendered abstract with the current APS/PRL contributor rule requiring no more than 600 characters including spaces.
- Cross-referenced every numbered Letter equation against the Supplemental citation.
- Inspected the decomposition code and recomputed raw and shifted closure residuals from both frozen decomposition CSV files.
- Visually inspected rendered pages containing Figs. 1 and 2 and checked the updated data-availability paragraph against bundle contents.

### Result

- The core numerical and scientific claims pass.
- The abstract-length test fails by roughly 47 characters.
- The Supplemental equation citation is unambiguously wrong.
- The archived implementation establishes the decomposition algebraically from independently contracted full/retained bilinear forms, but does not perform the stated direct `Hq` contraction; the quoted all-basis shifted closure maximum is also too small by a factor of about 2.6.

### Conclusion

- Status: contradicted — H1 is false for the exact final artifacts; H2 is supported.
- Confidence and evidence boundary: High confidence for the three local blockers and high confidence that they do not undermine the reported energies or the paper's central physical conclusion.
- Required next experiment: Shorten and recount the abstract; repair the Supplemental cross-reference; either add an independent tail contraction or rewrite the decomposition verification accurately and report a bound of at least `6.1e-14 Eh`; then recompile and repeat these three checks.

### Reflection

The successful numerical gates substantially strengthen the work, but they cannot validate journal-format text or prose that describes a stronger computational cross-check than the code actually performs. The final verdict is therefore a narrowly scoped revision, not a rejection of the science.

## Loop 003 — round-3 verification

### Observation

- The revised combined PDF has SHA-256 `7546e43371b20d4ad15568be3183cd4b96c18e39e6664f52aa2335ad87acc726`, ten letter-size pages, embedded fonts, and clean rendering of the revised Letter and Supplemental Material.
- Text extracted from the compiled PDF gives an abstract length of 548 characters including normalized spaces, below PRL's 600-character limit.
- Supplemental Eqs. `eq:s-rotation-h` and `eq:s-rotation-g` now separately define the one- and two-electron rotations at `supplement.tex:76-81`; `supplement.tex:243-245` cites exactly those equations.
- `main.tex:230-233` and `supplement.tex:372-379` now identify the decomposition as an algebraic reconstruction from directly contracted quantities and explicitly state that no separate `Hq` contraction was made.
- Recalculation from both frozen decomposition CSV files gives maximum normalized-energy recombination discrepancies of `6.039613253960852e-14 Eh` in the tracker basis and `2.275957200481571e-14 Eh` in the natural basis. The manuscript's rounded bound `6.1e-14 Eh` is conservative, and `figures/build_figures.py:318` enforces the same threshold for the plotted tracker-basis decomposition.

### Hypothesis

- H3: The three round-3 revisions remove every remaining substantive blocker without altering the frozen scientific results.
- Prediction if true: Each prior blocker passes its predeclared check, the compiled artifact remains structurally valid, and no new contradiction appears in the revised passages.
- Falsifier: Abstract length above 600, unresolved/misdirected equation references, a decomposition claim stronger than the implementation, a reported bound below the frozen maximum, or a broken PDF.

### Experiment

- Verified the exact PDF hash, extracted and normalized the rendered abstract, cross-referenced labels and citations in source, recomputed both-basis decomposition discrepancies from the frozen CSV files, inspected the figure-generation gate, checked page/font metadata and unresolved-reference text, and rendered representative revised PDF pages.

### Result

- All round-3 predictions passed. No quantitative result, figure, scope statement, or reproducibility claim was weakened or contradicted by the revisions.

### Conclusion

- Status: H3 proven within the stated artifact-review scope.
- Confidence and evidence boundary: High confidence for scientific, numerical, reproducibility, and PRL-format readiness of the exact reviewed anonymous artifact. Final author metadata and permanent public accession remain disclosed submission logistics rather than scientific blockers.
- Remaining uncertainty: None substantive within the authorized review scope.
- Referee verdict: `PUBLISHABLE`.

### Reflection

The revised wording now distinguishes direct Hamiltonian contractions from algebraic reconstruction precisely, while the conservative numerical threshold is tied to the actual frozen maximum. The prior objections are closed rather than merely rephrased.
