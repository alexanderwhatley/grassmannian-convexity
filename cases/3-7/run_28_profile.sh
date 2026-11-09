#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "$0")" && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf -- "$build_dir"' EXIT
cxx="${CXX:-c++}"
"$cxx" -O3 -std=c++17 -pthread "$root/move_refinement/verify_3x6_runs_independent.cpp" -o "$build_dir/verify_runs"
"$cxx" -O3 -std=c++17 -pthread "$root/move_refinement/verify_3x6_one_three_rank.cpp" -o "$build_dir/verify_rank"
python3 "$root/move_refinement/verify_3x6_separation.py"
"$build_dir/verify_runs"
cd "$root/move_refinement"
"$build_dir/verify_rank"
