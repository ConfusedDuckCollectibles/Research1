#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source_dir=$PWD
lean=$(printenv LEAN || echo lean)
lean_path=$(command -v "$lean")
version=$("$lean_path" --version)
[[ "$version" == *"version 4.19.0"* && "$version" == *"commit 6caaee842e94"* ]] ||
  { echo "Expected Lean 4.19.0 commit 6caaee842e94: $version" >&2; exit 1; }
printf '%s\n' "$version"

build=$(mktemp -d)
trap 'rm -rf "$build"' EXIT
export LEAN_PATH="$build"
export lean_path
jobs=$(printenv JOBS || echo 4)

c++ -O3 -std=c++17 -Wall -Wextra -Wpedantic "$source_dir/generate.cpp" -o "$build/gen"
"$build/gen" > "$build/certificates.json"
(cd "$build" && sha256sum --check "$source_dir/corpus.sha256")
SUBSET_CERTIFICATES="$build/certificates.json" python3 "$source_dir/audit_data.py"
python3 "$source_dir/emit_lean.py" --corpus "$build/certificates.json" --output "$build"
cp "$source_dir/Checker.lean" "$source_dir/Rows.lean" "$build/"
cd "$build"
"$lean_path" -o Checker.olean Checker.lean
"$lean_path" -o Rows.olean Rows.lean
"$lean_path" -o Templates.olean Templates.lean
"$lean_path" -o Controls.olean Controls.lean
echo 'Kernel controls passed (valid / truncated / extended / invalid-index certificates).'

# Every batch has an Eq.refl theorem for its cached run and a separate
# compiled cross-check against the pure-Lean row reconstruction.
xargs -P "$jobs" -n 1 bash -c \
  '"$lean_path" -o "$1.olean" "$1.lean" > "$1.log" 2>&1 || { cat "$1.log"; exit 1; }' -- \
  < modules.txt
while IFS= read -r batch; do cat "$batch.log"; done < modules.txt
"$lean_path" -o Manifest.olean Manifest.lean
python3 - <<'PY'
from pathlib import Path
mods=Path('modules.txt').read_text().splitlines()
assert len(mods)==81, len(mods)
counts={9:0,10:0,11:0}
for mod in mods:
    lines=Path(mod+'.log').read_text().splitlines()
    assert len(lines)==2 and lines[0]==f"'SubsetCover.{mod}.refuted' depends on axioms: [propext, Quot.sound]",lines
    word, name, count = lines[1].split()
    assert word=='AGREE' and name==mod
    counts[int(mod.split('_')[0][5:])]+=int(count)
assert counts=={9:62,10:3755,11:3755},counts
print('PASS:',len(mods),'kernel theorem batches;',sum(counts.values()),
      'concrete cached certificates kernel checked and separately matched by compiled Lean.')
PY
