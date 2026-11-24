#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"; src=$PWD
lean=$(command -v "$(printenv LEAN || echo lean)")
build=$(mktemp -d)
trap 'rm -rf "$build"' EXIT
export LEAN_PATH="$build"
c++ -O3 -std=c++17 -Wall -Wextra -Wpedantic "$src/generate.cpp" -o "$build/gen"
"$build/gen" > "$build/certificates.json"
(cd "$build" && sha256sum --check "$src/corpus.sha256")
python3 emit_lean.py --corpus "$build/certificates.json" --output "$build" --trial 1
python3 bench_agreement.py "$build/AgreementBench.lean"
cp Checker.lean Rows.lean "$build/"
cd "$build"
"$lean" -o Checker.olean Checker.lean
"$lean" -o Rows.olean Rows.lean
"$lean" -o Templates.olean Templates.lean
grep '^--' AgreementBench.lean
time timeout 70s "$lean" -o AgreementBench.olean AgreementBench.lean
