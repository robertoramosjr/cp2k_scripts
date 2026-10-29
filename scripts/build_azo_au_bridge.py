#!/usr/bin/env python3
"""Build azobenzene-between-two-Au-tips bridges (cis/trans x chemisorbed/physisorbed).

Electrodes: the two halves of a constricted Au nanowire (VASP POSCAR, split at its
single-atom neck); each half keeps its own internal geometry and is rigidly translated
so that its apex atom lands on the anchor point. The molecule comes from a TURBOMOLE
`coord` file (Bohr). Transport axis = z.

Anchoring:
  chem  apex Au binds the N lone pair (in the C-N=N plane, outward bisector) at --d-au-n.
        The tip-tip vector is then set by the two lone pairs and is aligned with z.
  phys  the centroid-to-centroid axis of the two rings is aligned with z; each apex sits
        --d-au-h from the outermost ring H along that axis (no covalent bond).

  anchored  4,4'-dissubstituted azo (--anchored LABEL=COORD, LABEL ends in _sh or _nh2). Anchor
        atoms = the two terminal S (thiolate: their H is REMOVED) or terminal N of -NH2 (H kept).
        Apex Au sits on the C->anchor extension at --d-au-s (S) or --d-au-n (N); the apex-apex
        vector is aligned with z. Assumptions: thiolate binding, on-top site.

Usage:
  build_azo_au_bridge.py --wire NanoWire_002_cart.vasp \
      --mol trans=azo_trans_hse06_m4_d3/coord --mol cis=azo_cis_hse06_m4_d3/coord \
      --outdir structures [--d-au-n 2.2] [--d-au-h 3.0] [--cell 22 22]
Writes <outdir>/<iso>_<anchor>.xyz (cell in the comment line) and prints a distance report.
"""
import argparse
from pathlib import Path

import numpy as np

BOHR = 0.529177210903


def read_poscar(path):
    lines = Path(path).read_text().split("\n")
    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[i].split()] for i in (2, 3, 4)]) * scale
    counts = [int(x) for x in lines[6].split()]
    k = 7
    if lines[k].strip().lower().startswith("s"):
        k += 1
    cart = lines[k].strip().lower().startswith(("c", "k"))
    k += 1
    xyz = np.array([[float(x) for x in lines[k + i].split()[:3]] for i in range(sum(counts))])
    return xyz if cart else xyz @ cell


def read_turbomole_coord(path):
    el, xyz = [], []
    for line in Path(path).read_text().split("\n")[1:]:
        if line.startswith("$"):
            break
        p = line.split()
        xyz.append([float(v) * BOHR for v in p[:3]])
        el.append(p[3].capitalize())
    return el, np.array(xyz)


def split_wire(au, tol=0.7):
    """Return (left, right) halves split at the central neck; apex = atom nearest the split."""
    au = au[np.argsort(au[:, 2])]
    layers = [[0]]
    for i in range(1, len(au)):
        if au[i, 2] - au[i - 1, 2] < tol:
            layers[-1].append(i)
        else:
            layers.append([i])
    half = len(layers) // 2
    left = au[[i for L in layers[:half] for i in L]]
    right = au[[i for L in layers[half:] for i in L]]
    return left, right, left[np.argmax(left[:, 2])], right[np.argmin(right[:, 2])]


def rotation_to_z(v):
    """Rotation matrix taking unit(v) to +z."""
    v = v / np.linalg.norm(v)
    z = np.array([0.0, 0.0, 1.0])
    axis = np.cross(v, z)
    s, c = np.linalg.norm(axis), float(v @ z)
    if s < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    a = axis / s
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def lone_pair_dirs(el, xyz):
    """Unit vectors of the two N lone pairs (outward bisector of the C-N=N angle)."""
    n_idx = [i for i, e in enumerate(el) if e == "N"]
    c_idx = [i for i, e in enumerate(el) if e == "C"]
    dirs = []
    for i in n_idx:
        j = [k for k in n_idx if k != i][0]
        c = min(c_idx, key=lambda k: np.linalg.norm(xyz[k] - xyz[i]))
        u = -((xyz[c] - xyz[i]) / np.linalg.norm(xyz[c] - xyz[i])
              + (xyz[j] - xyz[i]) / np.linalg.norm(xyz[j] - xyz[i]))
        dirs.append(u / np.linalg.norm(u))
    return n_idx, np.array(dirs)


def ring_centroids(el, xyz):
    """Centroids of the two phenyl rings (6 C each, split by nearest-N)."""
    n_idx = [i for i, e in enumerate(el) if e == "N"]
    c_idx = [i for i, e in enumerate(el) if e == "C"]
    groups = {n: [] for n in n_idx}
    for c in c_idx:
        groups[min(n_idx, key=lambda n: np.linalg.norm(xyz[c] - xyz[n]))].append(c)
    return [xyz[g].mean(axis=0) for g in groups.values()]


def anchor_points(el, xyz, mode, d_au_n, d_au_h):
    """Return (rotated molecule xyz, apex1, apex2) with the transport axis along z."""
    if mode == "chem":
        n_idx, u = lone_pair_dirs(el, xyz)
        t = [xyz[n] + d_au_n * ui for n, ui in zip(n_idx, u)]
        axis = t[1] - t[0]
    else:
        c1, c2 = ring_centroids(el, xyz)
        axis = c2 - c1
    rot = rotation_to_z(axis)
    xyz = (xyz - xyz.mean(axis=0)) @ rot.T
    if mode == "chem":
        n_idx, u = lone_pair_dirs(el, xyz)
        a1, a2 = (xyz[n] + d_au_n * ui for n, ui in zip(n_idx, u))
    else:
        h = [i for i, e in enumerate(el) if e == "H"]
        lo, hi = min(h, key=lambda i: xyz[i, 2]), max(h, key=lambda i: xyz[i, 2])
        a1 = xyz[lo] - np.array([0, 0, d_au_h])
        a2 = xyz[hi] + np.array([0, 0, d_au_h])
    if a1[2] > a2[2]:
        a1, a2 = a2, a1
    return xyz, a1, a2


def anchored_mol(el, xyz, grp, d_au_s, d_au_n):
    """Return (el, rotated xyz, apex1, apex2) for a -SH / -NH2 4,4'-dissubstituted azo."""
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    cov = {"H": 0.6, "C": 0.8, "N": 0.8, "S": 1.05}
    bonded = lambda i: [j for j in range(len(el)) if j != i and d[i, j] < 1.25 * (cov[el[i]] + cov[el[j]]) * 1.0 + 0.05]
    if grp == "sh":
        anc = [i for i, e in enumerate(el) if e == "S"]
        drop = [j for i in anc for j in bonded(i) if el[j] == "H"]
        dist = d_au_s
    else:
        anc = [i for i, e in enumerate(el) if e == "N" and sum(el[j] == "H" for j in bonded(i)) == 2]
        drop, dist = [], d_au_n
    assert len(anc) == 2, f"expected 2 anchor atoms, found {anc}"
    ups = []
    for i in anc:
        c = [j for j in bonded(i) if el[j] == "C"]
        assert len(c) == 1, f"anchor {i} bonded to C {c}"
        v = xyz[i] - xyz[c[0]]
        ups.append(v / np.linalg.norm(v))
    t = [xyz[i] + dist * u for i, u in zip(anc, ups)]
    rot = rotation_to_z(t[1] - t[0])
    cen = xyz.mean(axis=0)
    xyz = (xyz - cen) @ rot.T
    a = [(t[k] - cen) @ rot.T for k in (0, 1)]
    if a[0][2] > a[1][2]:
        a.reverse()
    keep = [i for i in range(len(el)) if i not in drop]
    return [el[i] for i in keep], xyz[keep], a[0], a[1]


def build(wire_au, el, xyz, mode, d_au_n, d_au_h):
    left, right, apex_l, apex_r = split_wire(wire_au)
    mol, a1, a2 = anchor_points(el, xyz, mode, d_au_n, d_au_h)
    left = left + (a1 - apex_l)
    right = right + (a2 - apex_r)
    return mol, left, right, a1, a2


def report(tag, el, mol, left, right, a1, a2):
    au = np.vstack([left, right])
    d = np.linalg.norm(au[:, None] - mol[None], axis=2)
    print(f"== {tag}: apex-apex {np.linalg.norm(a2 - a1):.2f} A | Au(left) {len(left)} Au(right) {len(right)}")
    for i, e in enumerate(el):
        if e == "N":
            print(f"   N{i}: nearest Au {d[:, i].min():.2f} A")
    heavy = [i for i, e in enumerate(el) if e != "H"]
    print(f"   min Au-C {d[:, [i for i in heavy if el[i] == 'C']].min():.2f} A | "
          f"min Au-H {d[:, [i for i, e in enumerate(el) if e == 'H']].min():.2f} A | "
          f"min Au-Au(left-right) {np.linalg.norm(left[:, None] - right[None], axis=2).min():.2f} A")


def write_xyz(path, el, mol, left, right, cell):
    allx = np.vstack([left, right, mol])
    allx[:, :2] += np.array(cell) / 2 - allx[:, :2].mean(axis=0)
    allx[:, 2] -= allx[:, 2].min() - 2.0
    lz = allx[:, 2].max() + 2.0
    syms = ["Au"] * (len(left) + len(right)) + el
    with open(path, "w") as f:
        f.write(f"{len(syms)}\nCell {cell[0]:.3f} {cell[1]:.3f} {lz:.3f} (provisional; Au {len(left)}+{len(right)}, mol {len(el)})\n")
        for s, p in zip(syms, allx):
            f.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wire", required=True, help="constricted Au nanowire POSCAR")
    ap.add_argument("--mol", action="append", default=[], metavar="ISO=COORD",
                    help="isomer label=TURBOMOLE coord (repeatable)")
    ap.add_argument("--anchored", action="append", default=[], metavar="LABEL=COORD",
                    help="e.g. cis_sh=anchored/cis_sh/coord (LABEL ends in _sh or _nh2)")
    ap.add_argument("--d-au-s", type=float, default=2.4, help="Au-S distance, anchored thiolate (A)")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--d-au-n", type=float, default=2.2, help="Au-N distance, chemisorbed (A)")
    ap.add_argument("--d-au-h", type=float, default=3.0, help="Au...H distance, physisorbed (A)")
    ap.add_argument("--cell", type=float, nargs=2, default=[22.0, 22.0], metavar=("LX", "LY"))
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    wire = read_poscar(a.wire)
    for spec in a.anchored:
        label, path = spec.split("=", 1)
        el, xyz = read_turbomole_coord(path)
        el, mol, a1, a2 = anchored_mol(el, xyz, label.rsplit("_", 1)[1], a.d_au_s, a.d_au_n)
        left, right, apex_l, apex_r = split_wire(wire)
        left, right = left + (a1 - apex_l), right + (a2 - apex_r)
        report(label, el, mol, left, right, a1, a2)
        dd = np.linalg.norm(np.vstack([left, right])[:, None] - mol[None], axis=2)
        for i, e in enumerate(el):
            if e in "SN" and e == ("S" if label.endswith("sh") else "N") and dd[:, i].min() < 3.0:
                print(f"   anchor {e}{i}: nearest Au {dd[:, i].min():.2f} A")
        write_xyz(a.outdir / f"{label}.xyz", el, mol, left, right, a.cell)
    for spec in a.mol:
        iso, path = spec.split("=", 1)
        el, xyz = read_turbomole_coord(path)
        for mode in ("chem", "phys"):
            mol, left, right, a1, a2 = build(wire, el, xyz, mode, a.d_au_n, a.d_au_h)
            tag = f"{iso}_{mode}"
            report(tag, el, mol, left, right, a1, a2)
            write_xyz(a.outdir / f"{tag}.xyz", el, mol, left, right, a.cell)


if __name__ == "__main__":
    main()
