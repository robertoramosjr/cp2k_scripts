#!/usr/bin/env python3
"""Build the single-anchor azobenzene half-bridge: one Au tip bonded to the molecule, the other
tip facing the free ring across a gap (plan v2 section 6, step A2).

The parent azobenzene (TURBOMOLE `coord`, Bohr) gets ONE para substituent (-SH thiolate or -NH2)
replacing the para H of ring A. Ring A sits at the bottom: C->anchor bond along -z, bottom apex Au
on its extension (--d-au-s / --d-au-n). The top apex is placed above the free para H of the TRANS
isomer at --d-gap; the same top-tip position is reused for the cis isomer (electrodes do not
move, only the isomer changes: that is the switch). Au tips are the two halves of the constricted
nanowire, rigid (make_junction_input.py --free-au 0 freezes all Au).

usage: build_azo_halfbridge.py --wire NanoWire_002_cart.vasp --anchor sh|nh2 \
           --trans <coord> --cis <coord> --outdir structures [--d-gap 3.0]
Writes <outdir>/{trans,cis}_<anchor>_half.xyz (cell in the comment line) and prints distances.
"""
import argparse
from pathlib import Path

import numpy as np

from build_azo_au_bridge import (read_poscar, read_turbomole_coord, rotation_to_z,
                                 split_wire, write_xyz)

COV = {"H": 0.31, "C": 0.76, "N": 0.71, "S": 1.05}


def neighbours(el, xyz, i, scale=1.25):
    return [j for j in range(len(el)) if j != i
            and np.linalg.norm(xyz[i] - xyz[j]) < scale * (COV[el[i]] + COV[el[j]])]


def rings(el, xyz):
    """Two phenyl rings as (ipso, para, ring atoms); ipso = ring C bonded to N, para = farthest ring C."""
    n_idx = [i for i, e in enumerate(el) if e == "N"]
    out = []
    for n in n_idx:
        ipso = [j for j in neighbours(el, xyz, n) if el[j] == "C"][0]
        ring, todo = {ipso}, [ipso]
        while todo:
            for j in neighbours(el, xyz, todo.pop()):
                if el[j] == "C" and j not in ring:
                    ring.add(j)
                    todo.append(j)
        para = max(ring, key=lambda j: np.linalg.norm(xyz[j] - xyz[ipso]))
        out.append((ipso, para, sorted(ring)))
    assert len(out) == 2 and all(len(r) == 6 for _, _, r in out), "ring detection failed"
    return out


def para_H(el, xyz, para):
    return [j for j in neighbours(el, xyz, para) if el[j] == "H"][0]


def substitute(el, xyz, anchor, ring_a):
    """Replace the para H of ring A by -S (thiolate, no H) or -NH2; return el, xyz, anchor index."""
    ipso, para, _ = ring_a
    h = para_H(el, xyz, para)
    u = (xyz[h] - xyz[para]) / np.linalg.norm(xyz[h] - xyz[para])
    keep = [i for i in range(len(el)) if i != h]
    el2, xyz2 = [el[i] for i in keep], xyz[keep]
    p2 = keep.index(para)
    if anchor == "sh":
        el2.append("S"); xyz2 = np.vstack([xyz2, xyz[para] + 1.77 * u])
        return el2, xyz2, len(el2) - 1, p2
    # -NH2: N in the ring plane, H's +-120/2 deg about the ring normal
    ring = xyz[ring_a[2]]
    nrm = np.linalg.svd(ring - ring.mean(0))[2][2]
    n_pos = xyz[para] + 1.38 * u
    w = np.cross(nrm, u)
    hs = [n_pos + 1.01 * (np.cos(np.radians(60)) * u + s * np.sin(np.radians(60)) * w) for s in (1, -1)]
    el2 += ["N", "H", "H"]
    xyz2 = np.vstack([xyz2, n_pos, hs[0], hs[1]])
    return el2, xyz2, len(el2) - 3, p2


def orient(el, xyz, anchor):
    """Substitute, then rotate/translate: anchor atom at the origin, its C->X bond along -z."""
    ra, rb = rings(el, xyz)
    el2, xyz2, a_idx, c_idx = substitute(el, xyz, anchor, ra)
    ra2, rb2 = rings(el2, xyz2)      # indices shifted by the removed H: redo on the new geometry
    if c_idx not in ra2[2]:          # ring A = the one holding the substituted para C
        ra2, rb2 = rb2, ra2
    v = xyz2[a_idx] - xyz2[c_idx]
    rot = rotation_to_z(-v)          # C->X onto -z
    xyz3 = (xyz2 - xyz2[a_idx]) @ rot.T
    h_free = para_H(el2, xyz3, rb2[1])
    return el2, xyz3, a_idx, h_free


def build(wire_au, anchor, el, xyz, d_anchor, d_gap, top_xy_z=None):
    el, mol, a_idx, h_free = orient(el, xyz, anchor)
    bottom = np.array([0.0, 0.0, -d_anchor])
    top = top_xy_z if top_xy_z is not None else mol[h_free] + np.array([0.0, 0.0, d_gap])
    left, right, apex_l, apex_r = split_wire(wire_au)
    left, right = left + (bottom - apex_l), right + (top - apex_r)
    return el, mol, left, right, top, h_free


def report(tag, el, mol, left, right, h_free):
    au = np.vstack([left, right])
    d = np.linalg.norm(au[:, None] - mol[None], axis=2)
    ser = lambda sym: d[:, [i for i, e in enumerate(el) if e == sym]].min()
    print(f"== {tag}: Au {len(left)}+{len(right)}, mol {len(el)} atoms, apex-apex {right[:, 2].min() - left[:, 2].max():.2f} A (z)")
    print(f"   free para H: z={mol[h_free, 2]:.2f}, nearest top Au {np.linalg.norm(right - mol[h_free], axis=1).min():.2f} A")
    print(f"   min Au-C {ser('C'):.2f} | Au-H {ser('H'):.2f} | top tip min to any mol atom "
          f"{np.linalg.norm(right[:, None] - mol[None], axis=2).min():.2f} A | "
          f"left-right Au-Au {np.linalg.norm(left[:, None] - right[None], axis=2).min():.2f} A")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wire", required=True)
    ap.add_argument("--anchor", choices=["sh", "nh2"], required=True)
    ap.add_argument("--trans", required=True, help="TURBOMOLE coord, trans azobenzene")
    ap.add_argument("--cis", required=True, help="TURBOMOLE coord, cis azobenzene")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--d-au-s", type=float, default=2.4, help="Au-S distance (A)")
    ap.add_argument("--d-au-n", type=float, default=2.2, help="Au-N distance (A)")
    ap.add_argument("--d-gap", type=float, default=3.0, help="top apex to free para H, trans (A)")
    ap.add_argument("--cell", type=float, nargs=2, default=[22.0, 22.0])
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    wire = read_poscar(a.wire)
    d_anchor = a.d_au_s if a.anchor == "sh" else a.d_au_n
    el_t, xyz_t = read_turbomole_coord(a.trans)
    el_c, xyz_c = read_turbomole_coord(a.cis)
    el, mol, left, right, top, h = build(wire, a.anchor, el_t, xyz_t, d_anchor, a.d_gap)
    report(f"trans_{a.anchor}", el, mol, left, right, h)
    write_xyz(a.outdir / f"trans_{a.anchor}_half.xyz", el, mol, left, right, a.cell)
    el, mol, left, right, _, h = build(wire, a.anchor, el_c, xyz_c, d_anchor, a.d_gap, top_xy_z=top)
    report(f"cis_{a.anchor}", el, mol, left, right, h)
    write_xyz(a.outdir / f"cis_{a.anchor}_half.xyz", el, mol, left, right, a.cell)


if __name__ == "__main__":
    main()
