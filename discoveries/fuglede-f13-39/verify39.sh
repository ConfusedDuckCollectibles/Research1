#!/bin/sh
# Exact primary replay; run with a C++17 compiler, Python 3.10+, and Node.js.
set -eu
cd "$(dirname "$0")"
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT HUP INT TERM
c++ -std=c++17 -O3 m3-planar-required.cpp -o "$d/planar"
c++ -std=c++17 -O3 m3-profile-grids.cpp -o "$d/profiles"
"$d/planar"
"$d/profiles" 6 9 > "$d/grids.out"
grep -E '^(L |TOTAL |doubles )' "$d/grids.out"
python3 m3-grid-orbits.py "$d/grids.out"
"$d/profiles" 7 7 | grep -E '^(L |TOTAL |doubles )'
python3 m3-small.py
python3 m3-skeleton.py
node finish39.js
