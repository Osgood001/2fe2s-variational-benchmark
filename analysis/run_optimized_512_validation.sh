#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
test "$(sha256sum 2fe_2s_30e_20o.fcidump | awk '{print $1}')" = "bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7"
test "$(sha256sum optimized_512_state.csv | awk '{print $1}')" = "1d9314c2c6e48fe82e45c932e8ca97f59d7156ec95ce7c19facad522f502a267"
mkdir -p output runtime-deps
python3 -m pip install --disable-pip-version-check --no-deps --target runtime-deps \
  "$h5py_wheel" "$wheel"
export PYTHONPATH="$PWD/runtime-deps${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64
python3 -u validate_optimized_512.py 2fe_2s_30e_20o.fcidump \
  optimized_512_state.csv reference_slater_condon.py \
  output/optimized_512_validation.json \
  2>&1 | tee output/validation.log
sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps
