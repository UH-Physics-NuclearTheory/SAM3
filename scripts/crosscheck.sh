#!/usr/bin/env bash
# Cross-check the C++ (cpp/sam3.h) and Mathematica (Mathematica/SAM-3.0.nb)
# SAM-3.0 implementations on identical random inputs.
#
#   EIGEN_INCLUDE=/path/to/eigen ./scripts/crosscheck.sh
#
# Requires: a C++17 compiler, Eigen, and wolframscript.
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"   # repo root (this script lives in scripts/)

EIGEN_INCLUDE="${EIGEN_INCLUDE:-}"
if [[ -z "$EIGEN_INCLUDE" ]]; then
  for c in /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3; do
    [[ -f "$c/Eigen/Dense" ]] && EIGEN_INCLUDE="$c" && break
  done
fi
[[ -f "$EIGEN_INCLUDE/Eigen/Dense" ]] || { echo "set EIGEN_INCLUDE to your Eigen include path"; exit 2; }
command -v wolframscript >/dev/null 2>&1 || { echo "wolframscript not found"; exit 2; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

echo "building bridge (sam3.h) ..."
g++ -std=c++17 -O2 "$HERE/cpp/bridge_sam3.cpp" -I"$HERE/cpp" -I"$EIGEN_INCLUDE" -o "$TMP/bridge_sam3"

# (d, s, Nmax): multi-charge (s=2,3) and higher order (3,4), plus multi-observable.
files=()
seed=1001
for cfg in "1 2 3" "1 3 4" "2 2 3" "1 3 3"; do
  read -r d s nm <<< "$cfg"
  out="$TMP/bridge_d${d}_s${s}_N${nm}.m"
  "$TMP/bridge_sam3" "$d" "$s" "$nm" "$seed" "$out"
  files+=("$out")
  seed=$((seed + 1))
done

echo "evaluating Mathematica KappaCE on the same inputs ..."
wolframscript -file "$HERE/Mathematica/crosscheck_cpp.wl" "${files[@]}"
