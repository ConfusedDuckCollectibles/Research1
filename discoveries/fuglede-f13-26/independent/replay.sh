#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
g++ -O3 -std=c++20 profile.cpp -o "$work/profiles"
g++ -O3 -std=c++20 caps.cpp -o "$work/caps"
{
  "$work/profiles"
  for testcase in '6 8' '7 7' '8 6' '8 5' '9 4'; do
    "$work/profiles" $testcase
  done
} | tee "$work/profile"
diff -u profile-replay.txt "$work/profile"
"$work/caps" $(awk '/^slope mask/ && $6=="[doubles" && $7==0 {print $3}' "$work/profile") | tee "$work/caps-out"
diff -u caps-replay.txt "$work/caps-out"
