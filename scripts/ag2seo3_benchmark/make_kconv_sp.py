#!/usr/bin/env python3
"""
make_kconv_sp.py -- k-mesh convergence by SINGLE POINTS (ENERGY_FORCE + analytical stress)
at ONE FIXED geometry, for every functional.

Why: the CELL_OPT k-mesh scan diverged for 6x3x5 and 8x4x6 (cell b grew past CELL_REF, SCF
garbage). Fixed-geometry single points isolate the k-sampling error from the optimizer:
dE/atom and dP (stress) vs the densest mesh. If the stress at the 5x3x4 minimum stays within
the CELL_OPT tolerance (100 bar) at denser meshes, 5x3x4 is converged.

The geometry is read from a CP2K .restart (top-level &CELL / &COORD; the nested &CELL_REF
is ignored). Same geometry for all functionals (k-error is nearly functional independent).

Usage:
  python make_kconv_sp.py --restart PBE/kmesh/k5_3_4/Ag2SeO3_PBE_k534-1.restart \
      --cutoff 500 --rel-cutoff 50 --out-dir kconv
"""
import argparse
import re
from pathlib import Path

from make_benchmark import FUNCTIONALS, block_input

MESHES = [(4, 2, 3), (5, 3, 4), (6, 3, 5), (8, 4, 6), (10, 5, 8)]


def read_restart(path):
    t = Path(path).read_text()
    i = t.index("&CELL\n", t.index("&SUBSYS"))
    cell = t[i:t.index("&END CELL\n", i)].split("&CELL_REF")[0]
    lat = [[float(x) for x in re.search(rf"^\s*{k}\s+(.*)$", cell, re.M).group(1).split()[:3]] for k in "ABC"]
    i = t.index("&COORD", t.index("&SUBSYS"))
    coord = t[i:t.index("&END COORD", i)]
    sp, fr = [], []
    for line in coord.splitlines()[1:]:
        w = line.split()
        if len(w) == 4 and w[0][0].isalpha() and w[0] != "SCALED":
            sp.append(w[0])
            fr.append([float(x) for x in w[1:]])
    assert sp and "SCALED T" in coord, "expected SCALED T coordinates in the restart"
    return lat, sp, fr


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--restart", required=True)
    p.add_argument("--cutoff", type=float, required=True)
    p.add_argument("--rel-cutoff", type=float, required=True)
    p.add_argument("--out-dir", default="kconv")
    p.add_argument("--functionals", nargs="+", default=list(FUNCTIONALS), choices=list(FUNCTIONALS))
    p.add_argument("--meshes", nargs="+", default=[".".join(map(str, m)) for m in MESHES],
                   help="e.g. 6.3.5 8.4.6")
    a = p.parse_args()
    lat, sp, fr = read_restart(a.restart)
    root = Path(a.out_dir)
    lst = []
    for f in a.functionals:
        for ms in a.meshes:
            m = tuple(int(x) for x in ms.split("."))
            d = root / f / ("k%d_%d_%d" % m)
            proj = "Ag2SeO3_%s_sp_k%d%d%d" % ((f,) + m)
            d.mkdir(parents=True, exist_ok=True)
            txt = block_input(proj, "ENERGY_FORCE", f, lat, sp, fr, a.cutoff, a.rel_cutoff, 0.01, m)
            txt = txt.replace("  METHOD Quickstep\n",
                              "  METHOD Quickstep\n  STRESS_TENSOR ANALYTICAL\n"
                              "  &PRINT\n    &STRESS_TENSOR ON\n    &END STRESS_TENSOR\n  &END PRINT\n", 1)
            (d / f"{proj}.inp").write_text(txt)
            lst.append(f"{d.resolve()} {proj}.inp")
    (root / "runs_kconv.txt").write_text("\n".join(lst) + "\n")
    print(f"{len(lst)} inputs; list: {root / 'runs_kconv.txt'}")


if __name__ == "__main__":
    main()
