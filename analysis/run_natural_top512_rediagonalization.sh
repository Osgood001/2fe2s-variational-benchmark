#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
test "$(sha256sum "$wheel" | awk '{print $1}')" = "37b0bccc55450311a55318cd643e851353331ddeab4fc0c0065e83c905e41502"
test "$(sha256sum "$h5py_wheel" | awk '{print $1}')" = "0cbd41f4e3761f150aa5b662df991868ca533872c95467216f2bec5fcad84882"
mkdir -p output runtime-deps
python3 -m pip install --disable-pip-version-check --no-deps --target runtime-deps \
  "$h5py_wheel" "$wheel"
export PYTHONPATH="$PWD/runtime-deps${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OPENBLAS_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64

python3 -u rediagonalize_topk_support.py \
  natural_orbital.fcidump natural_basis_ritz_vector.npy \
  reference_slater_condon.py output/natural_top512_support_rediagonalized.csv \
  output/natural_top512_rediagonalization.json \
  --k 512 \
  --expected-fcidump-sha256 adb401872023da5b0345ac473a121961d0553f9ddaaa808072423b0a6b24b2aa \
  --expected-vector-sha256 3855dcb7ce8c44b0f6d9c99fc5787aa2db339002b4d50ec2af8a814c45372d77 \
  --expected-reference-sha256 b655e63a9619983c6c7a30e3a82d47fb79bc35f93b7c41e1d6ae7772bc777f42 \
  --reference-topk-weight 0.7685072562843145 \
  --reference-fixed-energy-eh -116.43662614444148 \
  2>&1 | tee output/rediagonalization.log

sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps
