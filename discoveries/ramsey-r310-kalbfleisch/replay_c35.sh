#!/bin/sh
# No third-party graph data or large proof downloads are needed.
# Usage: ./replay_c35.sh /path/to/kissat /path/to/drat-trim
set -eu
[ "$#" -eq 2 ] || { printf '%s\n' "usage: $0 /path/to/kissat /path/to/drat-trim" >&2; exit 2; }
kissat=$1
trim=$2
python3 -B c35_extension.py --skip-exoo-baseline
python3 -B audit_c35_formulas.py --formulas-only
for b in pentagon bipartite star; do
    proof="$b.replayed.drat"
    rm -f "$proof"
    code=0
    "$kissat" --quiet --no-binary "$b.cnf" "$proof" >"$b.replayed.solver.log" 2>&1 || code=$?
    if [ "$code" -ne 20 ] || ! grep -q '^s UNSATISFIABLE$' "$b.replayed.solver.log"; then
        printf '%s\n' "UNRESOLVED: $b (exit $code)" >&2
        exit 1
    fi
    "$trim" "$b.cnf" "$proof" -q >"$b.replayed.verifier.log" 2>&1
    if ! tr '\r' '\n' <"$b.replayed.verifier.log" | grep -q '^s VERIFIED$'; then
        printf '%s\n' "UNVERIFIED: $b" >&2
        exit 1
    fi
    printf '%s\n' "VERIFIED: $b"
done
