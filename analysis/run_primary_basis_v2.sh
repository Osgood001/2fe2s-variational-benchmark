#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
fcidump_sha="bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7"
vector_sha="45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945"

test "$(sha256sum "$wheel" | awk '{print $1}')" = "37b0bccc55450311a55318cd643e851353331ddeab4fc0c0065e83c905e41502"
test "$(sha256sum "$h5py_wheel" | awk '{print $1}')" = "0cbd41f4e3761f150aa5b662df991868ca533872c95467216f2bec5fcad84882"
test "$(sha256sum 2fe_2s_30e_20o.fcidump | awk '{print $1}')" = "$fcidump_sha"
test "$(sha256sum fci_ritz_vector.npy | awk '{print $1}')" = "$vector_sha"

mkdir -p output runtime-deps
python3 -m pip install --disable-pip-version-check --no-deps --target runtime-deps \
  "$h5py_wheel" "$wheel"
export PYTHONPATH="$PWD/runtime-deps${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64

{
  date -u +'%Y-%m-%dT%H:%M:%SZ'
  uname -a
  python3 --version
  lscpu
  free -h
  sha256sum "$wheel" "$h5py_wheel" 2fe_2s_30e_20o.fcidump \
    fci_ritz_vector.npy analyze_state_compressibility.py run_primary_basis_v2.sh
} > output/environment.txt 2>&1

python3 -u analyze_state_compressibility.py \
  2fe_2s_30e_20o.fcidump fci_ritz_vector.npy output \
  --basis-label "tracker-orbital basis" \
  --expected-fcidump-sha256 "$fcidump_sha" \
  --expected-vector-sha256 "$vector_sha" \
  --reference-energy-eh -116.60560912042631 \
  --reference-energy-tolerance-eh 2e-10 \
  --max-full-residual-eh 1e-6 \
  --rhf-energy-eh -116.20581611838442 \
  --write-natural-orbital-fcidump \
  2>&1 | tee output/analysis.log

sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps
