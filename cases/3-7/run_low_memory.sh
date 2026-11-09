#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "$0")" && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf -- "$build_dir"' EXIT
mkdir -p "$build_dir/run" "$build_dir/grassmann_37_continuation/history_20261108"
gzip -dc "$root/low_memory/rank_by_sorted_key.bin.gz" > "$build_dir/grassmann_37_continuation/history_20261108/rank_by_sorted_key.bin"
gzip -dc "$root/low_memory/K_by_sorted_key.bin.gz" > "$build_dir/grassmann_37_continuation/history_20261108/K_by_sorted_key.bin"
"${CXX:-c++}" -O3 -std=c++17 -pthread "$root/low_memory/verify_low_memory.cpp" -o "$build_dir/verify_low_memory"
cd "$build_dir/run"
"$build_dir/verify_low_memory"
