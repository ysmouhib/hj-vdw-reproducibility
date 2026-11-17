#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mkdir -p "$ROOT/build"
CXX=${CXX:-g++}
for name in rabsearch_exhaustive rabung_weight_local simplex_warm_weighted vdw_tabu; do
    "$CXX" -O3 -std=c++17 "$ROOT/src/$name.cpp" -o "$ROOT/build/$name"
done
printf '%s\n' "Built search programs in $ROOT/build"
