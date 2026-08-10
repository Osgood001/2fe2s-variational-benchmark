# Hamiltonian-guided 512-determinant control

These files preserve the scientific construction code for the independent
512-determinant control used in Fig. 2. The reported energy is not trusted from
this code: `../validate_optimized_512.py` recomputes the frozen CSV with both an
explicit fermionic Slater--Condon implementation and a full-space PySCF
embedding.

The construction has two stages:

1. `exp_structure.py` grows an 80,000-determinant candidate space from the
   closed-shell determinant with Hamiltonian connections and perturbative
   addition estimates, then writes `cipsi.npz`.
2. `opt3.py` initializes a 512-determinant support from that pool, performs
   exact dense subspace diagonalization, ranks additions and removals with
   secular equations, accepts only energy-lowering swaps, and uses seeded
   large-neighborhood restarts.

Example commands for one deterministic seed and fixed wall-time budget are:

```bash
python exp_structure.py 2fe_2s_30e_20o.fcidump \
  --grow 80000 --budget 512 --outdir guided_pool

python opt3.py 2fe_2s_30e_20o.fcidump \
  --pool guided_pool/cipsi.npz --budget 512 --seed 11 \
  --ntry 64 --hours 3 --outdir guided_seed11
```

Because the stopping rule is a wall-time budget, the exact terminal support can
depend on hardware speed. The publication therefore freezes the selected state
as `../../data/optimized_512_state.csv` (SHA-256
`1d9314c2c6e48fe82e45c932e8ca97f59d7156ec95ce7c19facad522f502a267`).
The construction is a finite heuristic search and is not represented as a
global optimum. The CSV plus independent validator reproduce every numerical
claim made from this control.
