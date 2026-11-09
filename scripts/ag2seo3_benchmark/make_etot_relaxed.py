#!/usr/bin/env python3
"""make_etot_relaxed.py -- 'Etot' sigma series (GEO_OPT, FIXED cell = ISIF=2) starting from the
cell+ions of each functional's own finished production CELL_OPT (.restart).
Usage: make_etot_relaxed.py --restart R --func F --cutoff C --rel-cutoff RC --kpoints a b c --out-root DIR"""
import argparse, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_benchmark import block_input
p = argparse.ArgumentParser()
p.add_argument("--restart", required=True); p.add_argument("--func", required=True)
p.add_argument("--cutoff", type=float, required=True); p.add_argument("--rel-cutoff", type=float, required=True)
p.add_argument("--kpoints", type=int, nargs=3, required=True)
p.add_argument("--sigmas", type=float, nargs="+", default=[0.01, 0.1, 0.5, 1, 2, 3])
p.add_argument("--out-root", required=True); p.add_argument("--list", required=True)
a = p.parse_args()
txt = Path(a.restart).read_text()
cb = re.search(r"&CELL\n(.*?)&CELL_REF", txt, re.S).group(1)
lat = [[float(x) for x in re.search(rf"^\s*{k}\s+(.*)$", cb, re.M).group(1).split()[:3]] for k in "ABC"]
coord = re.search(r"&COORD\n(.*?)&END COORD", txt, re.S).group(1).splitlines()
assert any(l.split()[:2] == ["SCALED", "T"] for l in coord if l.strip())
sp, fr = [], []
for l in coord:
    t = l.split()
    if len(t) == 4: sp.append(t[0]); fr.append([float(x) for x in t[1:]])
lines = []
for s in a.sigmas:
    d = Path(a.out_root) / a.func / f"sigma{s:g}"; d.mkdir(parents=True, exist_ok=True)
    proj = f"Ag2SeO3_{a.func}_sigma{s:g}"
    (d / f"{proj}.inp").write_text(block_input(proj, "GEO_OPT", a.func, lat, sp, fr, a.cutoff, a.rel_cutoff, s, tuple(a.kpoints), False))
    lines.append(f"{d.resolve()} {proj}.inp")
with open(a.list, "a") as f: f.write("\n".join(lines) + "\n")
print(a.func, "natoms", len(sp), "|a|,|b|,|c| =", [round(sum(x*x for x in v)**.5, 4) for v in lat])
