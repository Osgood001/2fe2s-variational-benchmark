# Quantum Advantage Tracker submission #238

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

Paper, code, figures, and compact numerical results:
https://github.com/Osgood001/2fe2s-variational-benchmark

The calculation starts from the public CAS(30e,20o) FCIDUMP and a deterministic
two-entry CI guess. Matrix-free Davidson minimization is performed in the
complete fixed-`M_S=0`, `(N_alpha,N_beta)=(15,15)` determinant space of
dimension `240,374,016`. Explicit vectors are passed through four checkpointed
stages (`110/2/16/160` updates). No HCI or DMRG wavefunction is used.

The final normalized vector has SHA-256
`45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945`.
A fresh Hamiltonian action gives `E=-116.60560912042631 Eh` and
`||Hc-Ec||=6.901557706123275e-7 Eh`; an independent RDM contraction gives
`-116.60560912042571 Eh`. The solver flag is `converged=false`, so the energy
is reported only as an explicit variational upper bound. It is `184.12
microEh` below submission #187 and agrees with the rounded 2017 DMRG value
`-116.6056091 Eh`; no literature-minimum claim is made.

![Energy comparison and Davidson convergence](https://raw.githubusercontent.com/Osgood001/2fe2s-variational-benchmark/main/figures/figure1.png)

![Basis-dependent determinant compression](https://raw.githubusercontent.com/Osgood001/2fe2s-variational-benchmark/main/figures/figure2.png)

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
