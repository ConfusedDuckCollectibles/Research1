#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
d="$(mktemp -d)"
trap 'rm -rf "$d"' EXIT
c++ -O3 -std=c++17 -Wextra -Wall census.cpp -o "$d/census"
for n in {1..6}; do "$d/census" "$n" > "$d/k-$n.txt"; done
python3 -u check.py "$d"/k-{1,2,3,4,5,6}.txt
