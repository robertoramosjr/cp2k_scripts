#!/usr/bin/env python3
"""parse_kconv_sp.py -- tabulate the fixed-geometry k-mesh single points (kconv/).
dE (meV/atom) and pressure (bar, = STRESS| 1/3 Trace as printed) vs the densest finished mesh, per functional.
Reads only the last 'PROGRAM STARTED AT' segment of each .out (CP2K appends)."""
import glob, re, sys
from pathlib import Path
HA = 27211.386  # meV/hartree
NAT = 24
root = Path(sys.argv[1] if len(sys.argv) > 1 else "kconv")
for fn in sorted(p.name for p in root.iterdir() if p.is_dir()):
    rows = []
    for out in glob.glob(str(root / fn / "k*" / "*.out")):
        t = open(out).read(); t = t[t.rfind("PROGRAM STARTED AT"):]
        mesh = Path(out).parent.name[1:].replace("_", "x")
        nk = tuple(map(int, mesh.split("x")))
        ok = "PROGRAM ENDED" in t and "SCF run NOT converged" not in t
        e = re.findall(r"ENERGY\| Total FORCE_EVAL.*?\]\s+(-?[\d.]+)", t)
        p = re.findall(r"STRESS\| 1/3 Trace\s+(-?[\d.E+-]+)", t)
        rows.append((nk[0]*nk[1]*nk[2], mesh, ok, float(e[-1]) if e else None, float(p[-1]) if p else None))
    rows.sort()
    done = [r for r in rows if r[2] and r[3] is not None]
    if not done:
        print(f"{fn}: nothing finished"); continue
    eref, pref = done[-1][3], done[-1][4]
    print(f"\n{fn}  (ref = {done[-1][1]})\n  mesh      status  E [Ha]             dE [meV/at]  P [bar]   dP [bar]")
    for _, m, ok, e, p in rows:
        if ok and e is not None:
            print(f"  {m:8s}  OK      {e:.8f}  {(e-eref)*HA/NAT:10.3f}  {p:8.1f}  {p-pref:8.1f}")
        else:
            print(f"  {m:8s}  {'NC/run' :6s}")
