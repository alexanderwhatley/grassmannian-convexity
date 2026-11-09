#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "$0")" && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf -- "$build_dir"' EXIT
"${CXX:-c++}" -O3 -std=c++17 -pthread "$root/verify_upper_bound.cpp" -o "$build_dir/verify_upper_bound"
python3 "$root/verify_sharpness.py"
python3 "$root/verify_moment.py" 18204 "$root/simple_moment_12.json"
python3 "$root/verify_compact_witness.py" "$root/compact_witness.json"
python3 "$root/positive_cell/verify_constructed_12.py" "$root/simple_moment_117.json"
python3 "$root/verify_compact_witness.py" "$root/short_witness.json"
python3 "$root/verify_compact_witness.py" "$root/optimal_witness.json"
python3 "$root/verify_compact_witness.py" "$root/optimal_witness_integer.json"
python3 "$root/verify_compact_witness.py" "$root/simple_moment_35.json"
python3 "$root/verify_compact_witness.py" "$root/simple_moment_28.json"
python3 "$root/verify_compact_witness.py" "$root/simple_moment_29.json"
python3 "$root/positive_cell/verify_constructed_12.py" "$root/simple_moment_55.json"
python3 "$root/positive_cell/verify_constructed_12.py" "$root/simple_moment_41.json"
python3 "$root/positive_cell/verify_constructed_12.py" "$root/simple_moment_40.json"
python3 "$root/positive_cell/small_3_5_sign_graph.py"
python3 "$root/positive_cell/verify_minimum_28_independent.py"
python3 "$root/positive_cell/verify_optimal_p345.py"
python3 "$root/positive_cell/verify_3x5_group_rank_jump.py"
python3 "$root/positive_cell/verify_rank2_cyclic_jump.py"
"$root/run_28_profile.sh"
python3 "$root/positive_cell/exact_ift_construct_rational.py" \
  "$root/positive_cell/matrix_simple.json" "$root/optimal_witness.json" \
  "$build_dir/constructed.json"
python3 "$root/positive_cell/verify_constructed_12.py" "$build_dir/constructed.json"
cd "$build_dir"
./verify_upper_bound
