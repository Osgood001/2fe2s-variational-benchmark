#!/usr/bin/env bash
set -euo pipefail

wheel="wheels/pyscf-2.14.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
h5py_wheel="wheels/h5py-3.14.0-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
source_fcidump_sha="bd6ef31773f7b082171f73faa6c98924948e2c250b6fbec74504d64361c8c6c7"
source_vector_sha="45e63bb63cd10e953260b6f5dabb8a1dda5518fbf5959002c23888dacb599945"
natural_fcidump_sha="${NATURAL_FCIDUMP_SHA256:?set NATURAL_FCIDUMP_SHA256}"
rotation_sha="${NATURAL_ROTATION_SHA256:?set NATURAL_ROTATION_SHA256}"
natural_vector_sha="3855dcb7ce8c44b0f6d9c99fc5787aa2db339002b4d50ec2af8a814c45372d77"
reference_occupations_sha="2e9632ff525337eae6329eac5cd5e39b71bd08af85291a335033296b5a732b7d"

test "$(sha256sum natural_occupations_reference.csv | awk '{print $1}')" = "$reference_occupations_sha"
test "$(sha256sum "$wheel" | awk '{print $1}')" = "37b0bccc55450311a55318cd643e851353331ddeab4fc0c0065e83c905e41502"
test "$(sha256sum "$h5py_wheel" | awk '{print $1}')" = "0cbd41f4e3761f150aa5b662df991868ca533872c95467216f2bec5fcad84882"
test "$(sha256sum 2fe_2s_30e_20o.fcidump | awk '{print $1}')" = "$source_fcidump_sha"
test "$(sha256sum fci_ritz_vector.npy | awk '{print $1}')" = "$source_vector_sha"
test "$(sha256sum natural_orbital.fcidump | awk '{print $1}')" = "$natural_fcidump_sha"
test "$(sha256sum natural_orbital_rotation.csv | awk '{print $1}')" = "$rotation_sha"

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
    fci_ritz_vector.npy natural_orbital.fcidump natural_orbital_rotation.csv \
    transform_ci_natural_basis.py analyze_state_compressibility.py \
    run_natural_basis_transform.sh
} > output/environment.txt 2>&1

python3 -u transform_ci_natural_basis.py \
  2fe_2s_30e_20o.fcidump natural_orbital.fcidump fci_ritz_vector.npy \
  natural_orbital_rotation.csv output/natural_basis_ritz_vector.npy \
  output/orbital_transform.json \
  --expected-source-fcidump-sha256 "$source_fcidump_sha" \
  --expected-natural-fcidump-sha256 "$natural_fcidump_sha" \
  --expected-source-vector-sha256 "$source_vector_sha" \
  --expected-rotation-sha256 "$rotation_sha" \
  --reference-energy-eh -116.60560912042631 \
  --reference-residual-eh 6.901557706123275e-7 \
  2>&1 | tee output/transform.log

transformed_vector_sha="$(python3 -c 'import json; print(json.load(open("output/orbital_transform.json"))["output_vector_sha256"])')"
test "$transformed_vector_sha" = "$natural_vector_sha"
python3 -u analyze_state_compressibility.py \
  natural_orbital.fcidump output/natural_basis_ritz_vector.npy output \
  --basis-label "natural-orbital basis" \
  --expected-fcidump-sha256 "$natural_fcidump_sha" \
  --expected-vector-sha256 "$natural_vector_sha" \
  --reference-energy-eh -116.60560912042631 \
  --reference-energy-tolerance-eh 2e-9 \
  --max-full-residual-eh 1e-6 \
  --occupation-trace-abs-tolerance 2e-8 \
  --reference-occupations-csv natural_occupations_reference.csv \
  --expected-reference-occupations-sha256 "$reference_occupations_sha" \
  --reference-occupations-max-abs-tolerance 2e-8 \
  --natural-basis-rdm1-max-abs-tolerance 2e-8 \
  2>&1 | tee output/analysis.log

sha256sum output/* > output/output_sha256.txt
rm -rf runtime-deps
