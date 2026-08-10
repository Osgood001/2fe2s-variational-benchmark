#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
mkdir -p output runtime-deps
python3 -m pip install --disable-pip-version-check --no-deps --target runtime-deps \
  "$h5py_wheel" "$wheel"
export PYTHONPATH="$PWD/runtime-deps${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64

python3 -u rediagonalize_topk_support.py \
  2fe_2s_30e_20o.fcidump fci_ritz_vector.npy \
  reference_slater_condon.py output/top512_support_rediagonalized.csv \
  output/top512_rediagonalization.json \
  --k 512 \
  --expected-fcidump-sha256 bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7 \
  --expected-vector-sha256 45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945 \
  --expected-reference-sha256 b655e63a9619983c6c7a30e3a82d47fb79bc35f93b7c41e1d6ae7772bc777f42 \
  --reference-topk-weight 0.19134587182188278 \
  --reference-fixed-energy-eh -116.00034344510804 \
  2>&1 | tee output/rediagonalization.log

sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps
