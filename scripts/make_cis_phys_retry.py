#!/usr/bin/env python3
"""Two retries of the collapsed cis 'phys' junction (azo_au_bridge, 2026-10-27).

  restr : original cis_phys.xyz + a SHAKE distance constraint, centroid of ring 35-40 to the
          upper apex Au15 held at TARGET A (the free relaxation pulled that ring to 2.3 A).
  far   : no constraint; electrodes pulled apart and the molecule re-centred so that the
          closest molecule-Au distance is >= MIN_GAP A, then relaxed freely.

usage: make_cis_phys_retry.py <azo_au_bridge dir>
Writes <dir>/cis_phys_{restr,far}/{cis_restr,cis_far}.inp (+ the far xyz); inputs come from
make_junction_input.py (CUTOFF 600 / REL_CUTOFF 50 from the scan of this system).
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

TARGET = 3.3
MIN_GAP = 3.7
RING = [35, 36, 37, 38, 39, 40]
APEX = 15
HERE = Path(__file__).resolve().parent


def read(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0])
    cell = [float(x) for x in L[1].split()[1:4]]
    return cell, [(l.split()[0], np.array(l.split()[1:4], float)) for l in L[2:2 + n]]


def gap(atoms):
    mol = np.array([p for e, p in atoms if e != "Au"])
    au = np.array([p for e, p in atoms if e == "Au"])
    return np.linalg.norm(mol[:, None] - au[None], axis=2).min()


def far(cell, atoms):
    zc = np.mean([p[2] for e, p in atoms if e != "Au"])
    for d in np.arange(0.0, 6.0, 0.05):
        out = [(e, p + (np.array([0, 0, d]) if e == "Au" and p[2] > zc else 0)
                + (np.array([0, 0, d / 2]) if e != "Au" else 0)) for e, p in atoms]
        if gap(out) >= MIN_GAP:
            return [cell[0], cell[1], cell[2] + d], out, d
    raise SystemExit("MIN_GAP not reached")


def write_xyz(path, cell, atoms):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\nCell {cell[0]:.3f} {cell[1]:.3f} {cell[2]:.3f} (cis far, MIN_GAP {MIN_GAP})\n")
        for e, p in atoms:
            f.write(f"{e:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


def gen(xyz, proj, outdir):
    subprocess.run([sys.executable, str(HERE / "make_junction_input.py"), "--xyz", str(xyz), "--project", proj,
                    "--outdir", str(outdir), "--cutoff", "600", "--rel-cutoff", "50"], check=True)


root = Path(sys.argv[1])
cell, atoms = read(root / "structures" / "cis_phys.xyz")

gen(root / "structures" / "cis_phys.xyz", "cis_restr", root / "cis_phys_restr")
inp = root / "cis_phys_restr" / "cis_restr.inp"
txt = inp.read_text()
colvar = f"""    &COLVAR
      &DISTANCE
        &POINT
          ATOMS {' '.join(map(str, RING))}
        &END POINT
        &POINT
          ATOMS {APEX}
        &END POINT
        POINTS 1 2
      &END DISTANCE
    &END COLVAR
"""
coll = f"""    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR
      TARGET [angstrom] {TARGET}
    &END COLLECTIVE
"""
assert txt.count("  &CONSTRAINT\n") == 1 and txt.count("  &SUBSYS\n") == 1 and txt.count("  &END CONSTRAINT\n") == 1
txt = txt.replace("  &SUBSYS\n", "  &SUBSYS\n" + colvar, 1)
txt = txt.replace("  &END CONSTRAINT\n", coll + "  &END CONSTRAINT\n", 1)
inp.write_text(txt)

ncell, natoms, d = far(cell, atoms)
(root / "cis_phys_far").mkdir(exist_ok=True)
write_xyz(root / "cis_phys_far" / "cis_far_start.xyz", ncell, natoms)
gen(root / "cis_phys_far" / "cis_far_start.xyz", "cis_far", root / "cis_phys_far")
print(f"far: electrodes +{d:.2f} A (cell z {ncell[2]:.3f}), min gap {gap(natoms):.2f} A")
