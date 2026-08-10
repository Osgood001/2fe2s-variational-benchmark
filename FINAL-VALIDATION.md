# Final validation record

Date: 2026-08-09 (Asia/Shanghai)

## Scientific result

- Hamiltonian: CAS(30e,20o), fixed `(N_alpha,N_beta)=(15,15)`.
- FCIDUMP SHA-256: `bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7`.
- Final vector SHA-256: `45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945`.
- Direct Rayleigh energy: `-116.60560912042631 Eh`.
- Residual norm: `6.901557706123275e-7 Eh`.
- RDM/direct difference: `5.9686e-13 Eh`; `<S^2>=2.8569746519844883e-9`.
- Interpretation: an explicit variational upper bound in the stated sector;
  no two-sided residual error bar and no claim of improvement over rounded
  DMRG.

## Frozen Bohrium replay

- Bohr ID: `20529398`; terminal phase `completed`; exit code `0`.
- Natural FCIDUMP/vector/reference-occupation hashes matched their frozen
  values before scientific work.
- All 18 analysis checks passed. The compact eight-file output manifest also
  passed `sha256sum -c` after download.
- Natural-basis `K=4,194,304`: `W_K=0.9997758503591009`, energy error
  `0.000348775354297004 Eh`, full-space residual `0.027167832183283697 Eh`.
- Natural-basis top-512 same-support rediagonalization:
  `-116.48626027879892 Eh`, coefficient-relaxation gain
  `0.049634134357276594 Eh`, reference/PySCF difference
  `5.68e-14 Eh`, projected residual `1.68e-13 Eh`, full-space residual
  `0.311254416665089 Eh`.

## Document build

- Built on Bohrium with REVTeX 4.2 and BibTeX.
- Letter: 4 US-letter pages; Supplemental Material: 6 pages; combined: 10
  pages.
- No unresolved citations/references, fatal TeX errors, or overfull boxes.
- All PDF fonts are embedded; all pages and both standalone figures were
  rendered and visually inspected.
- The abstract extracted from the compiled Letter is 548 characters including
  spaces, below the PRL limit of 600.
- Letter SHA-256:
  `fd99101fc73d66c324c9f58273c6837fa34058020e171a13b72489fa67542bff`.
- Supplemental Material SHA-256:
  `c159d54988513bff4d5e2df4af82b23b41f977f257f4355e547fcf02942c7e2d`.
- Combined PDF SHA-256:
  `7546e43371b20d4ad15568be3183cd4b96c18e39e6664f52aa2335ad87acc726`.
- Figure 1 SHA-256:
  `05e66489caded6364ce01a430e7ad01ebdba3ed13b43afcb1d3f83893c1ec00c`.
- Figure 2 SHA-256:
  `ecb21cbfa030bc99e13bf77cf4be0644ec6ba423b528b8165ebc5e850068eed5`.

## Independent manuscript review

- The same referee reviewed every revision round.
- Exact final verdict: `PUBLISHABLE`.
- The reviewed combined artifact is the PDF with SHA-256
  `7546e43371b20d4ad15568be3183cd4b96c18e39e6664f52aa2335ad87acc726`.

## Deliberately pending human fields

- Approved author list and affiliations.
- Permission and durable public accession for the 1,922,992,256-byte final
  vector.
- Final public-incumbent recheck immediately before submitting the Tracker
  issue.

These administrative fields are not inferred or fabricated in the manuscript
or issue draft.
