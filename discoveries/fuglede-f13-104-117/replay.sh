#!/bin/sh
set -eu
cd "$(dirname "$0")"
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
python3 planar.py
c++ -std=c++17 -O2 -Wall -Wextra independent.cpp -o "$tmpdir/checker"
"$tmpdir/checker"
