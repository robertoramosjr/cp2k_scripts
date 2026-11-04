#!/usr/bin/env python3
"""
analyze_rdf.py -- crystallinity from RDF peaks for the Ag2SeO3 sigma scan.

For every <func>/etot/sigma<S>/ it reads the last frame of the GEO_OPT trajectory
(*-pos-1.xyz), builds the partial RDFs g_ab(r) with periodic images in the fixed
cell (POSCAR), and reports per partial: peak height g_max, peak position and FWHM
inside a window, plus the RMSD of the final positions against the initial POSCAR
(minimum image). A sharper/higher peak = more crystalline; broad/low = disordered.

All g(r) use the SAME bin width and cell, so heights are comparable across runs.
Windows (Angstrom) are the first-shell defaults and can be changed with --window.

Usage:
  python3 analyze_rdf.py --poscar ../vasp_inputs/POSCAR_Convencional --root . \
      --out rdf_summary.csv [--rmax 8 --dr 0.02 --sigma-smooth 0.0]
Writes rdf_summary.csv and rdf_curves/<func>_sigma<S>.dat (r, g_AgAg, g_AgSe, ...).
"""
import argparse
import csv
import re
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np

WINDOWS = {"Se-O": (1.4, 2.2), "Ag-O": (1.9, 3.0), "Ag-Se": (2.8, 4.2), "Ag-Ag": (2.6, 3.8),
           "O-O": (2.2, 3.4)}


def read_poscar(path):
    L = Path(path).read_text().splitlines()
    sc = float(L[1].split()[0])
    lat = np.array([[float(x) * sc for x in L[i].split()[:3]] for i in (2, 3, 4)])
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    n = sum(cnt)
    frac = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(n)])
    species = [e for e, c in zip(els, cnt) for _ in range(c)]
    return lat, species, frac


def last_frame(xyz):
    """Return (species, cart[N,3]) of the last frame of an XYZ trajectory."""
    L = Path(xyz).read_text().splitlines()
    n = int(L[0].split()[0])
    step = n + 2
    nfr = len(L) // step
    blk = L[(nfr - 1) * step + 2:(nfr - 1) * step + 2 + n]
    sp = [b.split()[0] for b in blk]
    xyz_ = np.array([[float(x) for x in b.split()[1:4]] for b in blk])
    return sp, xyz_, nfr


def images(lat, nmax=(2, 1, 2)):
    rng = [range(-m, m + 1) for m in nmax]
    return np.array([[i, j, k] for i in rng[0] for j in rng[1] for k in rng[2]]) @ lat


def rdf(cart, species, lat, a, b, rmax, dr, imgs):
    ia = [i for i, s in enumerate(species) if s == a]
    ib = [i for i, s in enumerate(species) if s == b]
    d = (cart[ib][None, :, None, :] + imgs[None, None, :, :] - cart[ia][:, None, None, :])
    d = np.linalg.norm(d, axis=-1).ravel()
    d = d[(d > 1e-6) & (d < rmax)]
    edges = np.arange(0, rmax + dr, dr)
    h, _ = np.histogram(d, edges)
    r = 0.5 * (edges[1:] + edges[:-1])
    vol = abs(np.linalg.det(lat))
    shell = 4 * np.pi * r ** 2 * dr
    rho_b = len(ib) / vol
    g = h / (len(ia) * shell * rho_b)
    return r, g


def peak(r, g, lo, hi):
    m = (r >= lo) & (r <= hi)
    if not m.any() or g[m].max() <= 0:
        return np.nan, np.nan, np.nan
    k = np.argmax(np.where(m, g, -1))
    half = g[k] / 2
    i = k
    while i > 0 and g[i] > half:
        i -= 1
    j = k
    while j < len(g) - 1 and g[j] > half:
        j += 1
    return g[k], r[k], r[j] - r[i]


def rmsd_min_image(cart, cart0, lat):
    inv = np.linalg.inv(lat)
    d = (cart - cart0) @ inv
    d -= np.round(d)
    return float(np.sqrt(np.mean(np.sum((d @ lat) ** 2, axis=1))))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--poscar", required=True)
    p.add_argument("--root", default=".")
    p.add_argument("--out", default="rdf_summary.csv")
    p.add_argument("--rmax", type=float, default=8.0)
    p.add_argument("--dr", type=float, default=0.02)
    a = p.parse_args()

    lat, sp0, frac = read_poscar(a.poscar)
    cart0 = frac @ lat
    imgs = images(lat)
    pairs = [f"{x}-{y}" for x, y in combinations_with_replacement(sorted(set(sp0)), 2)]
    rows = []
    cdir = Path(a.root) / "rdf_curves"
    cdir.mkdir(exist_ok=True)
    for xyz in sorted(Path(a.root).glob("*/etot/sigma*/*-pos-1.xyz")):
        func = xyz.parts[-4]
        sig = float(re.search(r"sigma([0-9.]+)", xyz.parts[-2]).group(1))
        out = next(xyz.parent.glob("*.out"), None)
        done = bool(out and "OPTIMIZATION COMPLETED" in out.read_text(errors="ignore"))
        sp, cart, nfr = last_frame(xyz)
        if sp != sp0:
            print(f"[skip] {xyz}: atom order differs from POSCAR")
            continue
        row = dict(func=func, sigma_eV=sig, frames=nfr, converged=done,
                   rmsd_A=round(rmsd_min_image(cart, cart0, lat), 4))
        curves = {}
        for pr in pairs:
            x, y = pr.split("-")
            r, g = rdf(cart, sp, lat, x, y, a.rmax, a.dr, imgs)
            curves[pr] = g
            lo, hi = WINDOWS.get(pr, (1.5, 4.0))
            h, rp, w = peak(r, g, lo, hi)
            row[f"{pr}_gmax"], row[f"{pr}_rpeak"], row[f"{pr}_fwhm"] = (round(h, 3), round(rp, 3), round(w, 3))
        np.savetxt(cdir / f"{func}_sigma{sig:g}.dat",
                   np.column_stack([r] + [curves[k] for k in pairs]), header="r " + " ".join(pairs))
        rows.append(row)
    rows.sort(key=lambda d: (d["func"], d["sigma_eV"]))
    if not rows:
        print("no trajectories found")
        return
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} runs -> {a.out}")


if __name__ == "__main__":
    main()
