#!/usr/bin/env python3
"""Emit three direct, non-cached concrete Q11 row-bucket equality theorems."""
import json
import sys
from pathlib import Path
from emit_lean import array

corpus=json.loads(Path(__file__).with_name("certificates.json").read_text())
template=next(t for t in corpus if t["dimension"]==11)
longest=max(template["entries"],key=lambda e:len(e[1]))
pool=[(ds, cert, s, rows) for ds,cert,cache in template["entries"] for s,rows in cache]
cases=[(longest[0],longest[1],s,rs)
       for s,rs in sorted(longest[2], key=lambda p:(-len(p[1]),p[0]))[:2]]
cases.append(max((c for c in pool if c[0]!=longest[0]),key=lambda c:len(c[3])))
lines=[
    "import Templates",
    "set_option maxRecDepth 8000",
    "set_option maxHeartbeats 0",
    "namespace SubsetCover.AgreementBench",
]
for i,(d,c,s,rs) in enumerate(cases):
    lines.append(f"-- {len(c)} preorder nodes; {len(rs)} rows in this site bucket.")
    data="["+(", ".join(map(array,rs)))+"]"
    lines.append(f"theorem bucket{i} : (t11.board {array(d)}).rowsAt {s} = {data} := by rfl")
    lines.append(f"#print axioms bucket{i}")
lines.append("end SubsetCover.AgreementBench")
Path(sys.argv[1]).write_text("\n".join(lines)+"\n")
