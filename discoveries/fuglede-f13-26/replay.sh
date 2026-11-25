#!/bin/sh
set -eu
cd "$(dirname "$0")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT INT TERM
c++ -O3 -std=c++17 profiles.cpp -o "$work/profiles"
{
    "$work/profiles" 6 8
    "$work/profiles" 7 7
    "$work/profiles" 8 6
    "$work/profiles" 8 5
    "$work/profiles" 9 4
} > "$work/result"
diff -u profiles.out "$work/result"
python3 controls.py "$work/result"
