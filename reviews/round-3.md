# Same-referee review round 3

Verdict: revision required.

## Referee blockers

1. The rendered abstract was 646 characters, above the PRL 600-character
   limit.
2. The natural-orbital section of the Supplemental Material cited Letter
   Eq. (2), which is the Rayleigh quotient rather than the integral rotation.
3. The retained--tail paragraph described a directly verified closure at
   `2.3e-14 Eh`, although the tail term was reconstructed by subtraction. A
   two-basis recomputation from the frozen files gave a maximum floating-point
   discrepancy of `6.0396132539608516e-14 Eh`.

## Revision

1. The abstract was shortened and recompiled. Text extracted from the final
   PDF is 548 characters.
2. The one- and two-electron rotations now carry separate Supplemental labels,
   and the natural-orbital paragraph cites those equations.
3. The Letter and Supplemental Material now state explicitly that the
   coupling and tail terms are algebraically reconstructed from directly
   contracted `E_Ritz`, `E_K`, and `<c|H|phi_K>`; no independent `Hq`
   contraction is implied. The reported recombination bound is rounded up to
   `6.1e-14 Eh`. The figure-generation gate was tightened to the same bound.

The revised Letter and Supplemental Material compile on Bohrium without
unresolved references, fatal errors, or overfull boxes and remain 4 and 6
pages, respectively.
