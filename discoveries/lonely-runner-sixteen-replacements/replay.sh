#!/bin/sh
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
"${CXX:-g++}" -std=c++17 -O3 -Wall -Wextra "$here/generate.cpp" -o "$tmp/generate"
case "$*" in
    "") set -- 2 3 4 ;;
    5) ;;
    *) echo "usage: sh replay.sh [5]" >&2; exit 2 ;;
esac
for k in "$@"; do
    "$tmp/generate" "$k" 18000000 "$tmp/certificate.txt"
    "${PYTHON:-python3}" "$here/verify.py" "$tmp/certificate.txt"
    rm "$tmp/certificate.txt"
done
