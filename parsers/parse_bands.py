#!/usr/bin/env python3
"""
parse_bands.py -- step 04: CP2K .bs file -> bands.csv, gap, plot.

Reads the band-structure file written by &PRINT%BAND_STRUCTURE:
  # Set 1: ...
  #  Special point 1  x y z  LABEL
  #  Point 1  Spin 1:  kx ky kz  weight
  #   Band    Energy [eV]     Occupation
         1    -41.12345       2.00000
k-point coordinates are fractional (B_VECTOR) and are converted to a
Cartesian path length with the reciprocal lattice of the cell in the .inp.
Energies are shifted so the VBM (highest state with occupation > --occ-tol)
is zero. Format checked on a real cp2k/2026.1 .bs (c-ZrO2, 2026-09-24):
one "# Set" per KPOINT_SET, special points evenly spaced by NPOINTS.

  python ~/work_cp2k/parsers/parse_bands.py --run-dir 04_bands
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_blocks as cb  # noqa: E402

RE_POINT = re.compile(r"#\s*Point\s+(\d+)\s+Spin\s+(\d+):\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)")
RE_SPECIAL = re.compile(r"#\s*Special point\s+\d+\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+(\S+)")
RE_SET = re.compile(r"#\s*Set\s+(\d+)")


def parse_bs(path: Path):
    sets, cur, kp = [], None, None
    for line in path.read_text().splitlines():
        if RE_SET.match(line):
            cur = dict(special=[], points=[])
            sets.append(cur)
        elif (m := RE_SPECIAL.match(line)) and cur is not None:
            cur["special"].append(([float(x) for x in m.groups()[:3]], m.group(4)))
        elif (m := RE_POINT.match(line)):
            kp = dict(spin=int(m.group(2)), k=[float(x) for x in m.groups()[2:5]], bands=[])
            cur["points"].append(kp)
        elif line.strip() and not line.lstrip().startswith("#") and kp is not None:
            parts = line.split()
            kp["bands"].append((float(parts[1]), float(parts[2])))
    return sets


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", default=".")
    p.add_argument("--occ-tol", type=float, default=1e-3)
    p.add_argument("--emin", type=float, default=-6.0)
    p.add_argument("--emax", type=float, default=8.0)
    a = p.parse_args()

    d = Path(a.run_dir)
    bs = sorted(d.glob("*.bs"))
    inp = sorted(d.glob("*.inp"))
    if not bs or not inp:
        sys.exit("ERROR: need <project>.bs and <project>.inp in --run-dir.")
    recip = cb.read_restart_structure(inp[0]).lattice.reciprocal_lattice.matrix
    sets = parse_bs(bs[0])
    xs, E, O, ticks, x0 = [], [], [], [], 0.0
    for s in sets:
        pts = [q for q in s["points"] if q["spin"] == 1]
        prev = None
        for q in pts:
            k = np.array(q["k"]) @ recip
            if prev is not None:
                x0 += np.linalg.norm(k - prev)
            prev = k
            xs.append(x0)
            E.append([b[0] for b in q["bands"]])
            O.append([b[1] for b in q["bands"]])
        nsp = len(s["special"])
        if pts and nsp > 1:  # special points are evenly spaced in index (NPOINTS per segment)
            step = (len(pts) - 1) / (nsp - 1)
            ticks += [(xs[-len(pts) + round(j * step)], s["special"][j][1]) for j in range(nsp)]
    E, O = np.array(E), np.array(O)
    vbm = E[O > a.occ_tol].max()
    cbm = E[O <= a.occ_tol].min() if (O <= a.occ_tol).any() else float("nan")
    ik_v = np.unravel_index(np.where(O > a.occ_tol, E, -np.inf).argmax(), E.shape)[0]
    ik_c = np.unravel_index(np.where(O <= a.occ_tol, E, np.inf).argmin(), E.shape)[0]
    kind = "direct" if ik_v == ik_c else "indirect"
    print(f"VBM {vbm:.4f} eV, CBM {cbm:.4f} eV, gap {cbm - vbm:.4f} eV ({kind})")
    with open(d / "bands.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x"] + [f"band{i + 1}" for i in range(E.shape[1])])
        for x, row in zip(xs, E - vbm):
            w.writerow([f"{x:.6f}"] + [f"{e:.5f}" for e in row])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(xs, E - vbm, c="k", lw=0.8)
        for x, lab in ticks:
            ax.axvline(x, c="gray", lw=0.5)
        ax.set_xticks([x for x, _ in ticks], [lab.replace("GAMMA", "Γ") for _, lab in ticks])
        ax.set_ylim(a.emin, a.emax)
        ax.set_xlim(xs[0], xs[-1])
        ax.set_ylabel("E - E_VBM [eV]")
        fig.tight_layout()
        fig.savefig(d / "bands.png", dpi=150)
    except Exception as exc:
        print(f"[INFO] plot skipped: {exc}")


if __name__ == "__main__":
    main()
