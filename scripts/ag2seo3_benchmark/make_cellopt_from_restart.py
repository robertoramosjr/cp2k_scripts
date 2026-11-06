#!/usr/bin/env python3
"""Build a CELL_OPT input from the cell+coords of a CP2K .restart, with CELL_REF = 1.15x the NEW cell.

Fixes the case where the cell grew beyond the original CELL_REF (fixed-grid reference) during a
previous CELL_OPT.  Reuses block_input() from make_benchmark.py.
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_benchmark import block_input  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("--restart", required=True)
p.add_argument("--func", required=True, choices=["PBE", "PBEsol", "revPBE", "BLYP"])
p.add_argument("--cutoff", type=float, required=True)
p.add_argument("--rel-cutoff", type=float, required=True)
p.add_argument("--kpoints", type=int, nargs=3, required=True)
p.add_argument("--sigma", type=float, default=0.01)
p.add_argument("--out-dir", required=True)
p.add_argument("--project", required=True)
a = p.parse_args()

txt = Path(a.restart).read_text()
cell_blk = re.search(r"&CELL\n(.*?)&CELL_REF", txt, re.S).group(1)  # only the real cell (before CELL_REF)
lat = [[float(x) for x in re.search(rf"^\s*{k}\s+(.*)$", cell_blk, re.M).group(1).split()[:3]] for k in "ABC"]
coord = re.search(r"&COORD\n(.*?)&END COORD", txt, re.S).group(1).splitlines()
assert any(l.split()[:2] == ["SCALED", "T"] for l in coord if l.strip()), "expected SCALED T in restart"
sp, fr = [], []
for l in coord:
    t = l.split()
    if len(t) == 4:
        sp.append(t[0]); fr.append([float(x) for x in t[1:]])
out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
(out / f"{a.project}.inp").write_text(block_input(
    a.project, "CELL_OPT", a.func, lat, sp, fr, a.cutoff, a.rel_cutoff, a.sigma, tuple(a.kpoints)))
print(f"{out.resolve()} {a.project}.inp   natoms={len(sp)} |b|={sum(x*x for x in lat[1])**.5:.3f}")
