#!/usr/bin/env python3
"""make_bands_sp.py -- single point (ENERGY) on a finished GEO_OPT geometry + &BAND_STRUCTURE eigenvalues along a
high-symmetry path of the cell's OWN reciprocal axes (P2_1/c, unique axis b: points in B_VECTOR of the input cell).
SCF mesh/sigma/cutoffs identical to the GEO_OPT, so occupations are the ones the relaxation saw.
Usage: make_bands_sp.py --restart R --func F --sigma S --cutoff C --rel-cutoff RC --kpoints a b c
                        --out-root DIR --list FILE [--npoints N]"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_benchmark import block_input
import re
P = {"GAMMA": (0, 0, 0), "X": (.5, 0, 0), "Y": (0, .5, 0), "Z": (0, 0, .5), "A": (.5, .5, 0),
     "B": (.5, 0, .5), "C": (0, .5, .5), "D": (.5, .5, .5)}
PATH = [["GAMMA", "X", "A", "Y", "GAMMA", "Z", "B", "D", "C", "Z"], ["X", "B"], ["Y", "C"], ["A", "D"]]
p = argparse.ArgumentParser()
p.add_argument("--restart", required=True); p.add_argument("--func", required=True)
p.add_argument("--sigma", type=float, required=True)
p.add_argument("--cutoff", type=float, required=True); p.add_argument("--rel-cutoff", type=float, required=True)
p.add_argument("--kpoints", type=int, nargs=3, required=True)
p.add_argument("--npoints", type=int, default=10)
p.add_argument("--out-root", required=True); p.add_argument("--list", required=True)
a = p.parse_args()
txt = Path(a.restart).read_text()
cb = re.search(r"&CELL\n(.*?)(?:&CELL_REF|&END CELL\n)", txt, re.S).group(1)
lat = [[float(x) for x in re.search(rf"^\s*{k}\s+(?:\[\w+\]\s+)?(.*)$", cb, re.M).group(1).split()[:3]] for k in "ABC"]
coord = re.search(r"&COORD\n(.*?)&END COORD", txt, re.S).group(1).splitlines()
assert any(l.split()[:2] == ["SCALED", "T"] for l in coord if l.strip())
sp, fr = [], []
for l in coord:
    t = l.split()
    if len(t) == 4: sp.append(t[0]); fr.append([float(x) for x in t[1:]])
proj = f"Ag2SeO3_{a.func}_sigma{a.sigma:g}_bands"
d = Path(a.out_root) / a.func / f"sigma{a.sigma:g}"; d.mkdir(parents=True, exist_ok=True)
inp = block_input(proj, "ENERGY", a.func, lat, sp, fr, a.cutoff, a.rel_cutoff, a.sigma, tuple(a.kpoints), False)
nmos = 208  # same number of bands printed for every sigma (SCF ADDED_MOS stays as in the GEO_OPT)
sets = ""
for br in PATH:
    sets += "        &KPOINT_SET\n          UNITS B_VECTOR\n"
    for n in br: sets += f"          SPECIAL_POINT {n} {P[n][0]:.8f} {P[n][1]:.8f} {P[n][2]:.8f}\n"
    sets += f"          NPOINTS {a.npoints}\n        &END KPOINT_SET\n"
dft_print = (f"    &PRINT\n      &BAND_STRUCTURE\n        FILE_NAME {proj}.bs\n        ADDED_MOS {nmos}\n"
             + sets + "      &END BAND_STRUCTURE\n    &END PRINT\n")
assert inp.count("  &END DFT\n") == 1 and "  &DFT\n" in inp
inp = inp.replace("  &END DFT\n", dft_print + "  &END DFT\n")
(d / f"{proj}.inp").write_text(inp)
with open(a.list, "a") as f: f.write(f"{d.resolve()} {proj}.inp\n")
print(a.func, a.sigma, "ADDED_MOS", nmos, "->", d / f"{proj}.inp")
