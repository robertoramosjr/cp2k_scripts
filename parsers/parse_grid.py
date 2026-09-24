#!/usr/bin/env python3
"""
parse_grid.py -- step 01: energy and Gaussian-per-grid distribution per point.

Criterion (CP2K how-to, 2D method):
  1. energy: |E(v) - E(v_max)| < --threshold-mev-atom for v and every larger v;
  2. grid:   the finest grid (level 1) is used (count > 0) while the MAJORITY
             of Gaussians sits on the coarser grids (sum levels >= 2 > level 1).
             Most Gaussians on grid 1 means REL_CUTOFF/CUTOFF too low for the
             multigrid to do its job; an empty grid 1 means wasted CUTOFF.
The recommended value is the smallest one passing both.

  python ~/work_cp2k/parsers/parse_grid.py --scan-dir 01_grid_convergence/cutoff
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_output as co  # noqa: E402
from core.cp2k_blocks import HARTREE_EV  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan-dir", required=True)
    p.add_argument("--threshold-mev-atom", type=float, default=1.0)
    p.add_argument("--plot", default="grid_convergence.png")
    a = p.parse_args()

    root = Path(a.scan_dir)
    meta = json.loads((root / "scan_meta.json").read_text())
    phase, nat = meta["phase"], meta["natoms"]
    rows = []
    for d in sorted((x for x in root.iterdir() if x.is_dir()),
                    key=lambda x: float(x.name.replace(phase, ""))):
        out = co.find_output(d)
        if not out:
            print(f"[WARNING] {d.name}: no output, skipped")
            continue
        t = co.read(out)
        e = co.energies_ha(t)
        if not e or not co.ended(t):
            print(f"[WARNING] {d.name}: no energy / did not end, skipped")
            continue
        g = co.grid_counts(t)
        rows.append(dict(value=float(d.name.replace(phase, "")), E_ha=e[-1], grids=g,
                         ranks=co.mpi_ranks(t)))
    if len(rows) < 2:
        sys.exit("ERROR: fewer than 2 finished points.")
    eref = rows[-1]["E_ha"]
    for r in rows:
        r["dE_mev_atom"] = abs(r["E_ha"] - eref) * HARTREE_EV * 1000 / nat
        c = [n for _, n, _ in r["grids"]]
        r["grid_ok"] = bool(c) and c[0] > 0 and sum(c[1:]) > c[0]
    rec = None
    for i, r in enumerate(rows):
        tail_ok = all(x["dE_mev_atom"] < a.threshold_mev_atom for x in rows[i:])
        if tail_ok and r["grid_ok"]:
            rec = r["value"]
            break
    print(f"{phase:>10} {'E [Ha]':>18} {'dE meV/at':>10}  counts per grid (1=finest)  ok")
    for r in rows:
        print(f"{r['value']:>10g} {r['E_ha']:>18.10f} {r['dE_mev_atom']:>10.3f}  "
              f"{[n for _, n, _ in r['grids']]}  {'grid' if r['grid_ok'] else '-'}")
    print(f"\nRecommended {phase.upper()}: {rec if rec is not None else 'NONE (extend the scan)'} "
          f"(threshold {a.threshold_mev_atom} meV/atom, reference = {rows[-1]['value']:g})")
    (root / "grid_summary.json").write_text(json.dumps(dict(phase=phase, recommended=rec, rows=rows,
                                                             threshold_mev_atom=a.threshold_mev_atom), indent=1))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.semilogy([r["value"] for r in rows[:-1]], [max(r["dE_mev_atom"], 1e-4) for r in rows[:-1]], "o-")
        ax.axhline(a.threshold_mev_atom, ls="--", c="gray")
        ax.set_xlabel(f"{phase.upper()} [Ry]")
        ax.set_ylabel("|E - E_ref| [meV/atom]")
        fig.tight_layout()
        fig.savefig(root / a.plot, dpi=150)
    except Exception as exc:  # plotting is optional
        print(f"[INFO] plot skipped: {exc}")


if __name__ == "__main__":
    main()
