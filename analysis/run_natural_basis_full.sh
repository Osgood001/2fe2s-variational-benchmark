#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
natural_fcidump_sha="${NATURAL_FCIDUMP_SHA256:?set NATURAL_FCIDUMP_SHA256 in the job command}"

test "$(sha256sum "$wheel" | awk '{print $1}')" = "37b0bccc55450311a55318cd643e851353331ddeab4fc0c0065e83c905e41502"
test "$(sha256sum "$h5py_wheel" | awk '{print $1}')" = "0cbd41f4e3761f150aa5b662df991868ca533872c95467216f2bec5fcad84882"
test "$(sha256sum natural_orbital.fcidump | awk '{print $1}')" = "$natural_fcidump_sha"

mkdir -p output runtime-deps work/A work/B1 work/B2 work/C
python3 -m pip install --disable-pip-version-check --no-deps --target runtime-deps \
  "$h5py_wheel" "$wheel"
export PYTHONPATH="$PWD/runtime-deps${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64
export PYSCF_MAX_MEMORY=220000
export EXPECTED_FCIDUMP_SHA256="$natural_fcidump_sha"

{
  date -u +'%Y-%m-%dT%H:%M:%SZ'
  uname -a
  python3 --version
  lscpu
  free -h
  sha256sum "$wheel" "$h5py_wheel" natural_orbital.fcidump \
    natural_orbital_rotation.csv initialize_fci_ritz.py \
    continue_fci_checkpointed.py analyze_state_compressibility.py \
    run_natural_basis_full.sh
} > output/environment.txt 2>&1

FCI_MAX_CYCLE=110 FCI_MAX_SPACE=14 \
  python3 -u initialize_fci_ritz.py natural_orbital.fcidump work/A \
  2>&1 | tee output/stage-A.log

FCI_MAX_CYCLE=2 FCI_MAX_SPACE=64 FCI_ENERGY_TOL=1e-13 \
FCI_RESIDUAL_TOL=1e-8 FCI_CHECKPOINT_EVERY=1 \
  python3 -u continue_fci_checkpointed.py natural_orbital.fcidump \
  work/A/fci_ritz_vector.npy work/B1 2>&1 | tee output/stage-B1.log
rm -f work/A/fci_ritz_vector.npy

FCI_MAX_CYCLE=16 FCI_MAX_SPACE=64 FCI_ENERGY_TOL=1e-13 \
FCI_RESIDUAL_TOL=1e-8 FCI_CHECKPOINT_EVERY=8 \
  python3 -u continue_fci_checkpointed.py natural_orbital.fcidump \
  work/B1/fci_ritz_vector.npy work/B2 2>&1 | tee output/stage-B2.log
rm -f work/B1/fci_ritz_vector.npy work/B1/latest_checkpoint.npy

FCI_MAX_CYCLE=160 FCI_MAX_SPACE=64 FCI_ENERGY_TOL=1e-13 \
FCI_RESIDUAL_TOL=1e-8 FCI_CHECKPOINT_EVERY=8 \
  python3 -u continue_fci_checkpointed.py natural_orbital.fcidump \
  work/B2/latest_checkpoint.npy work/C 2>&1 | tee output/stage-C.log
rm -f work/B2/fci_ritz_vector.npy work/B2/latest_checkpoint.npy

mv work/C/fci_ritz_vector.npy output/natural_basis_ritz_vector.npy
cp work/C/continuation_result.json output/natural_basis_continuation_result.json
vector_sha="$(sha256sum output/natural_basis_ritz_vector.npy | awk '{print $1}')"
python3 -u analyze_state_compressibility.py \
  natural_orbital.fcidump output/natural_basis_ritz_vector.npy output \
  --basis-label "natural-orbital basis" \
  --expected-fcidump-sha256 "$natural_fcidump_sha" \
  --expected-vector-sha256 "$vector_sha" \
  --reference-energy-eh -116.60560912042631 \
  --reference-energy-tolerance-eh 2e-6 \
  --max-full-residual-eh 2e-6 \
  2>&1 | tee output/analysis.log

rm -f work/C/latest_checkpoint.npy
sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps work
