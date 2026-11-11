#!/usr/bin/env python3
"""make_cellopt_exp.py -- production CELL_OPT (ISIF=3 equivalent) per functional, starting from the EXPERIMENTAL
Ag2SeO3 structure (ICSD 78388 = COD 1509694). Reuses the grid/k-mesh already converged for the same material
(CUTOFF 500, REL_CUTOFF 50; k 5x3x4 for PBE/PBEsol/revPBE, 8x4x6 for BLYP). Sigma 0.01 eV.
CELL_REF = 1.15x the initial cell, except 1.5x for BLYP/revPBE (their cells swelled 22-50 % in earlier runs).
Usage: make_cellopt_exp.py --poscar P --out-root DIR --list FILE"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_benchmark import block_input, read_poscar
CFG = {"PBE": ((5, 3, 4), 1.15), "PBEsol": ((5, 3, 4), 1.15), "revPBE": ((5, 3, 4), 1.5), "BLYP": ((8, 4, 6), 1.5)}
p = argparse.ArgumentParser()
p.add_argument("--poscar", required=True); p.add_argument("--out-root", required=True); p.add_argument("--list", required=True)
p.add_argument("--cutoff", type=float, default=500); p.add_argument("--rel-cutoff", type=float, default=50)
a = p.parse_args()
lat, sp, fr = read_poscar(a.poscar)
lines = []
for f, (k, cr) in CFG.items():
    d = Path(a.out_root) / f / "prod"; d.mkdir(parents=True, exist_ok=True)
    proj = f"Ag2SeO3_{f}_exp_prod"
    (d / f"{proj}.inp").write_text(block_input(proj, "CELL_OPT", f, lat, sp, fr, a.cutoff, a.rel_cutoff, 0.01, k, False, cell_ref_factor=cr))
    lines.append(f"{d.resolve()} {proj}.inp")
Path(a.list).write_text("\n".join(lines) + "\n"); print("\n".join(lines))
