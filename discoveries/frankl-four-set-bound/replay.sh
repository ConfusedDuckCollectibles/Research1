#!/bin/sh
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
CXX=${CXX:-c++}
python3 - "$here/liu-catalog8.txt" <<'PY'
from hashlib import sha256
from pathlib import Path
import sys
p=Path(sys.argv[1])
if sha256(p.read_bytes()).hexdigest()!="15dac91a167ccd9482d406f7cf054f9b6facdc23a826764618a4d4db865ffcd1":
    raise SystemExit("Unexpected bytes in Liu's pinned catalogue")
PY
"$CXX" -O3 -std=c++20 "$here/cards15.cpp" -o "$tmp/cards"
"$CXX" -O3 -std=c++20 "$here/verify24.cpp" -o "$tmp/verify24"
"$tmp/cards" "$here/liu-catalog8.txt" "$tmp/catalog.txt" > "$tmp/cards.txt"
python3 "$here/prepare24.py" "$tmp/cards.txt" "$tmp/catalog.txt" "$tmp/links.txt"
"$tmp/verify24" "$tmp/catalog.txt" "$tmp/links.txt"
printf '%s\n' 'PASS: FC(4,10) <= 24 from the stated upstream inputs.'
