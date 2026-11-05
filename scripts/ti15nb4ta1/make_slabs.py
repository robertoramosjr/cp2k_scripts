#!/usr/bin/env python3
"""
make_slabs.py -- slab geometries from the 20-atom Nb4Ta1Ti15 BCC supercell (numpy only).

Miller indices refer to the CUBIC BCC axes. The supercell of the CONTCAR is
  A = -2x, B = 2y - z, C = -y - 2z  (cubic axes, units of a)   det = -10 -> 10 cubic cells
so the cubic frame is recovered as E = inv(S0) @ L (rows = cubic vectors in Cartesian, A).
The alloy is disordered: the periodicity is the 20-atom lattice L, NOT the BCC lattice, hence
 * in-plane cell = primitive 2D sublattice of L inside the plane (hkl),
 * one "period" along the normal = smallest L-translation leaving the plane family (spacing d),
   always holding exactly 20 atoms,
 * terminations = cuts through the largest height gaps inside one period.
Slab = N periods, c = (thickness + vacuum) along the normal (orthogonal cell).

  python3 make_slabs.py --bulk bulk/04_cellopt/bulk_relaxed.vasp --planes "1,1,0;1,0,1;..." --nperiods 2 3 4
"""
import argparse, itertools, math
from pathlib import Path
import numpy as np

S0 = np.array([[-2, 0, 0], [0, 2, -1], [0, -1, -2]])


def read_poscar(path):
    L = Path(path).read_text().splitlines()
    sc = float(L[1].split()[0])
    lat = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * sc
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    n = sum(cnt)
    frac = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(n)])
    sp = [e for e, c in zip(els, cnt) for _ in range(c)]
    return lat, sp, frac


def write_poscar(path, lat, sp, cart, comment):
    order = []
    for e in sp:
        if e not in order:
            order.append(e)
    inv = np.linalg.inv(lat)
    with open(path, "w") as f:
        f.write(comment + "\n1.0\n")
        for v in lat:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write(" ".join(order) + "\n" + " ".join(str(sp.count(e)) for e in order) + "\nDirect\n")
        for e in order:
            for s, r in zip(sp, cart):
                if s == e:
                    f.write("  %.10f %.10f %.10f\n" % tuple(r @ inv))


def plane_cell(L, E, hkl):
    """u1, u2 (in-plane, integer combos of L rows), c3 (with n.c3 = g), all in Cartesian."""
    n = np.array(hkl)
    p = S0 @ n                      # n . (m S0) = m . p
    g = math.gcd(*[int(abs(x)) for x in p]) if any(p) else 1
    rng = range(-7, 8)
    ker = [np.array(m) for m in itertools.product(rng, rng, rng) if any(m) and np.dot(p, m) == 0]
    cart = lambda m: m @ L
    ker.sort(key=lambda m: np.linalg.norm(cart(m)))
    G = np.linalg.inv(E).T @ n      # reciprocal (plane normal), 1/A
    area_expect = abs(np.linalg.det(L)) * np.linalg.norm(G) / g
    u = None
    for i, m1 in enumerate(ker[:60]):
        for m2 in ker[i + 1:80]:
            a = np.linalg.norm(np.cross(cart(m1), cart(m2)))
            if abs(a - area_expect) < 1e-6 * area_expect + 1e-6:
                u = (m1, m2)
                break
        if u:
            break
    assert u, f"no primitive 2D cell found for {hkl}"
    m3 = next(np.array(m) for m in sorted(itertools.product(rng, rng, rng),
                                          key=lambda m: np.linalg.norm(cart(np.array(m))))
              if np.dot(p, m) == g)
    return cart(u[0]), cart(u[1]), cart(m3), g


def build(L, sp, fr, hkl, E, nper, term, vac):
    u1, u2, c3, g = plane_cell(L, E, hkl)
    nh = np.cross(u1, u2)
    nh = nh / np.linalg.norm(nh)
    if nh @ c3 < 0:
        nh = -nh
    d = nh @ c3                                  # period along the normal
    P = np.array([u1, u2, nh])                   # in-plane frame (orthonormalising the 3rd axis)
    inv2 = np.linalg.inv(np.array([u1, u2, nh]))
    r0 = fr @ L
    h = r0 @ nh
    h0 = h.min()
    r0 = r0 - np.floor((h - h0) / d - 1e-12)[:, None] * c3   # fold every atom into one period [h0, h0+d)
    h = r0 @ nh
    # gaps inside one period
    hm = (h - h0) % d
    order = np.argsort(hm)
    hs = hm[order]
    gaps = np.diff(np.append(hs, hs[0] + d))
    gap_idx = np.argsort(-gaps)[:max(term, 1)]
    k = gap_idx[term - 1] if term - 1 < len(gap_idx) else gap_idx[-1]
    z0 = hs[(k + 1) % len(hs)] if k + 1 < len(hs) else hs[0] + d
    z0 = z0 - 1e-3                                # cut just below the first atom after the gap
    atoms, spp = [], []
    for t in range(-1, nper + 1):
        for r, s in zip(r0, sp):
            rr = r + t * c3
            hh = rr @ nh - h0 - z0
            if 0 <= hh < nper * d - 1e-9:
                U = np.array([u1, u2])
                frac2 = np.linalg.solve(U @ U.T, U @ rr)
                frac2 -= np.floor(frac2)
                pos = frac2[0] * u1 + frac2[1] * u2 + hh * nh
                atoms.append(pos)
                spp.append(s)
    atoms = np.array(atoms)
    assert len(atoms) == 20 * nper, (hkl, len(atoms), nper)
    thick = nper * d
    cell = np.array([u1, u2, (thick + vac) * nh])
    atoms = atoms + 0.5 * vac * nh                 # centre the slab in the vacuum
    return cell, spp, atoms, d, np.linalg.norm(np.cross(u1, u2)), gaps[np.argsort(-gaps)][:term]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bulk", default="VASP_references/CONTCAR (1)")
    ap.add_argument("--planes", default="1,1,0;1,0,1;0,1,1;1,-1,0;1,0,-1;0,1,-1;1,0,0;0,1,0;0,0,1;1,1,1;1,1,-1;1,-1,1;-1,1,1")
    ap.add_argument("--nperiods", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--terminations", type=int, default=1, help="number of terminations per plane (largest gaps)")
    ap.add_argument("--vacuum", type=float, default=25.0)
    ap.add_argument("--out", default="slabs")
    a = ap.parse_args()
    L, sp, fr = read_poscar(a.bulk)
    E = np.linalg.inv(S0) @ L
    print("cubic axes lengths [A]:", np.round(np.linalg.norm(E, axis=1), 4), " right-handed:", np.linalg.det(E) > 0)
    print(f"{'hkl':>8} {'g':>2} {'d [A]':>6} {'area [A2]':>9} {'atoms/period':>12}  " + "  ".join(f"N={n}: thk/atoms" for n in a.nperiods))
    for hs in a.planes.split(";"):
        hkl = tuple(int(x) for x in hs.split(","))
        row = []
        for n in a.nperiods:
            for t in range(1, a.terminations + 1):
                cell, spp, at, d, area, gp = build(L, sp, fr, hkl, E, n, t, a.vacuum)
                tag = "".join(str(x).replace("-", "m") for x in hkl)
                out = Path(a.out) / f"hkl_{tag}"
                out.mkdir(parents=True, exist_ok=True)
                write_poscar(out / f"slab_{tag}_N{n}_t{t}.vasp", cell, spp, at, f"Nb4Ta1Ti15 slab ({hs}) N={n} term={t}")
            row.append(f"{n * d:5.1f} Å/{len(spp):3d}")
        print(f"{hs:>8} {math.gcd(*[int(abs(x)) for x in S0 @ np.array(hkl)]):>2} {d:6.2f} {area:9.1f} {20:>12}  " + "  ".join(f"{r:>14}" for r in row))


if __name__ == "__main__":
    main()
