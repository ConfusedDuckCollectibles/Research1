#!/usr/bin/env python3
"""Emit kernel-checkable theorem modules; no third-party Python packages."""
import argparse
import json
from pathlib import Path

def array(xs):
    return "[" + ", ".join(map(str, xs)) + "]"

def entry(e, *, pure):
    fields=[array(e[0]), array(e[1])]
    if not pure:
        fields.append("[" + ", ".join(f"({site}, [{', '.join(map(array, rows))}])"
                       for site,rows in e[2]) + "]")
    return "⟨" + ", ".join(fields) + "⟩"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--corpus", type=Path, default=Path("certificates.json"))
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--batch", type=int, default=96)
    p.add_argument("--trial", type=int)
    p.add_argument("--pure", action="store_true")
    a=p.parse_args()
    assert a.batch>0
    corpus=json.loads(a.corpus.read_text())
    dst=a.output
    dst.mkdir(parents=True,exist_ok=True)
    intro="import Rows\n\nnamespace SubsetCover\n"
    configs=intro
    modules=[]
    controls=[]
    for c in corpus:
        dim, h, w, trim, k, bans, entries=[c[key] for key in
            ("dimension","h","width","trim","tail","bans","entries")]
        assert len(bans)==h+w+1
        configs += f"def t{dim} : Template := ⟨{h}, {w}, {trim}, {k}, #{array(bans)}⟩\n"
        e=entries[0]
        controls.append(f"def ctl{dim} : Entry := {entry(e,pure=True)}")
        controls.append(f"example : passed t{dim} ctl{dim} = true := by rfl")
        controls.append(f"example : passed t{dim} {{ctl{dim} with certificate := {array(e[1][:-1])}}} = false := by rfl")
        controls.append(f"example : passed t{dim} {{ctl{dim} with certificate := {array(e[1]+[0])}}} = false := by rfl")
        controls.append(f"example : passed t{dim} {{ctl{dim} with certificate := {array([h+w+1]+e[1][1:])}}} = false := by rfl")
        if a.trial is not None:
            entries=entries[:a.trial]
        for num,start in enumerate(range(0,len(entries),a.batch)):
            name=f"Batch{dim}_{num:03}"
            modules.append(name)
            batch=entries[start:start+a.batch]
            es=",\n  ".join(entry(e,pure=a.pure) for e in batch)
            typ="Entry" if a.pure else "CachedEntry"
            pred="passed" if a.pure else "acceptedCached"
            board=(f"(t{dim}.board e.deficits)" if a.pure else f"(e.board t{dim})")
            sound="finite_batch_sound" if a.pure else "cached_batch_sound"
            compiled=""
            if not a.pure:
                compiled=(f"#eval (show IO Unit from do\n"
                          f"  unless entries.all (fun e => e.agrees t{dim}) do\n"
                          f"    throw (IO.userError \"cached rows differ from independent Lean reconstruction\")\n"
                          f"  IO.println \"AGREE {name} {len(batch)}\")\n")
            contents=("import Templates\n\n"
                "set_option maxRecDepth 8000\n"
                "set_option maxHeartbeats 0\n\n"
                f"namespace SubsetCover.{name}\n"
                f"def entries : List {typ} := [\n  {es}\n]\n\n"
                f"theorem accepted : entries.all ({pred} t{dim}) = true := by rfl\n"
                f"theorem refuted (e : {typ}) (h : e ∈ entries) : ¬ Residual {board} [] :=\n"
                f"  {sound} accepted h\n"
                "#print axioms refuted\n"
                + compiled +
                f"end SubsetCover.{name}\n")
            (dst/f"{name}.lean").write_text(contents)
    (dst/"Templates.lean").write_text(configs+"end SubsetCover\n")
    (dst/"Controls.lean").write_text(
        "import Templates\n\nset_option maxRecDepth 8000\nset_option maxHeartbeats 0\n\nnamespace SubsetCover\n"
        + "\n".join(controls) + "\n\nend SubsetCover\n")
    (dst/"Manifest.lean").write_text("".join(f"import {m}\n" for m in modules)
        + "import Controls\n\n#print axioms SubsetCover.check_sound\n")
    (dst/"modules.txt").write_text("\n".join(modules)+"\n")

if __name__=="__main__":
    main()
