#!/usr/bin/env python3
"""analyze_bands_occ.py -- gap, Gaussian-smeared occupations and conduction electrons from the SP band runs.
The .bs 'Occupation' column is integer filling (Ne = 184.000 exactly), NOT the SCF smearing, so occupations are
rebuilt as f = 0.5*erfc((e - EF)/sigma) per spin with EF read from the .out (Hartree -> eV). Check: N(EF) over the
path must be ~NELEC. Usage: analyze_bands_occ.py ROOT [--csv out.csv] [--nocc 92]   (ROOT/<F>/sigma<S>/*.bs,*.out)"""
import argparse, csv, glob, math, re
from pathlib import Path
HA = 27.211386245988
p = argparse.ArgumentParser(); p.add_argument("root"); p.add_argument("--csv"); p.add_argument("--nocc", type=int, default=92)
a = p.parse_args()
def bands(bs):
    ks, cur = [], None
    for l in open(bs):
        if l.startswith("#  Point"): cur = []; ks.append(cur)
        elif l.startswith("#") or not l.strip(): continue
        else: t = l.split(); cur.append((int(t[0]), float(t[1])))
    return ks
rows = []
for fd in sorted(Path(a.root).iterdir()):
    if not fd.is_dir(): continue
    for sd in sorted(fd.glob("sigma*"), key=lambda x: float(x.name[5:])):
        bs, out = glob.glob(f"{sd}/*.bs"), glob.glob(f"{sd}/*.out")
        if not bs or not out: continue
        sg = float(sd.name[5:]); ks = bands(bs[0])
        seg = open(out[0], errors="ignore").read().split("PROGRAM STARTED AT")[-1]
        ef = float(re.findall(r"Fermi energy:\s+(-?\d+\.\d+)", seg)[-1]) * HA
        f = lambda e: 0.5 * math.erfc((e - ef) / sg)
        vb = [max(e for i, e in k if i <= a.nocc) for k in ks]; cb = [min(e for i, e in k if i > a.nocc) for k in ks]
        g0 = [e for i, e in ks[0] if i == a.nocc + 1][0]
        rows.append(dict(func=fd.name, sigma=sg, nk=len(ks), nb=len(ks[0]), ef=ef, vbm=max(vb), cbm=min(cb),
                         gap=min(cb) - max(vb), gap_direct=min(c - v for v, c in zip(vb, cb)),
                         n_path=sum(2 * f(e) for k in ks for i, e in k) / len(ks),
                         f_cb1_gamma=f(g0), n_cond=sum(2 * f(e) for k in ks for i, e in k if i > a.nocc) / len(ks)))
hdr = "func  sigma   EF     VBM    CBM    gap   gap_dir  N(path)  f(b93,G)  e-_cond"
print(hdr)
for r in rows:
    print(f"{r['func']:6s}{r['sigma']:6g} {r['ef']:6.2f} {r['vbm']:6.2f} {r['cbm']:6.2f} {r['gap']:6.2f} {r['gap_direct']:6.2f}"
          f"  {r['n_path']:8.3f}  {r['f_cb1_gamma']:.3f}   {r['n_cond']:6.2f}")
if a.csv:
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
