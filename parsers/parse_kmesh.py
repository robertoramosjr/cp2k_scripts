#!/usr/bin/env python3
"""
parse_kmesh.py -- step 02: double criterion on relaxed energy and volume.

Only meshes whose CELL_OPT printed OPTIMIZATION COMPLETED are used. The
recommended mesh is the sparsest one for which it AND every denser mesh are
within --threshold-mev-atom (energy per atom) and --threshold-vol-pct
(relaxed volume) of the densest converged mesh. Also prints its k-density
n_i*|a_i| -- the number steps 03/07 reuse on other cells.

  python ~/work_cp2k/parsers/parse_kmesh.py --scan-dir 02_kmesh_convergence
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from core import cp2k_blocks as cb  # noqa: E402
from core import cp2k_output as co  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan-dir", required=True)
    p.add_argument("--threshold-mev-atom", type=float, default=1.0)
    p.add_argument("--threshold-vol-pct", type=float, default=0.5)
    a = p.parse_args()

    root = Path(a.scan_dir)
    meta = json.loads((root / "scan_meta.json").read_text())
    nat = meta["natoms"]
    rows = []
    for d in sorted(x for x in root.iterdir() if x.is_dir() and x.name.startswith("k")):
        mesh = [int(x) for x in d.name[1:].split("-")]
        out = co.find_output(d)
        t = co.read(out) if out else ""
        if not co.opt_converged(t):
            print(f"[WARNING] {d.name}: CELL_OPT not converged, skipped")
            continue
        rst = co.find_restart(d)
        st = cb.read_structure(rst) if rst else None
        rows.append(dict(mesh=mesh, nk=mesh[0] * mesh[1] * mesh[2], E_ev=co.final_energy_ev(t),
                         V=st.volume if st else co.last_volume(t),
                         density=cb.kmesh_density(st.lattice, mesh) if st else None))
    if len(rows) < 2:
        sys.exit("ERROR: fewer than 2 converged meshes.")
    rows.sort(key=lambda r: r["nk"])
    ref = rows[-1]
    for r in rows:
        r["dE_mev_atom"] = abs(r["E_ev"] - ref["E_ev"]) * 1000 / nat
        r["dV_pct"] = abs(r["V"] - ref["V"]) / ref["V"] * 100
    rec = None
    for i, r in enumerate(rows):
        if all(x["dE_mev_atom"] < a.threshold_mev_atom and x["dV_pct"] < a.threshold_vol_pct for x in rows[i:]):
            rec = r
            break
    print(f"{'mesh':>10} {'E [eV]':>16} {'dE meV/at':>10} {'V [A3]':>10} {'dV %':>7} {'dens':>6}")
    for r in rows:
        print(f"{'x'.join(map(str, r['mesh'])):>10} {r['E_ev']:>16.6f} {r['dE_mev_atom']:>10.3f} "
              f"{r['V']:>10.3f} {r['dV_pct']:>7.3f} {r['density'] or 0:>6.1f}")
    if rec:
        print(f"\nRecommended mesh {'x'.join(map(str, rec['mesh']))}, k-density {rec['density']:.1f} A")
    else:
        print("\nNo mesh passes both criteria -- add denser meshes.")
    (root / "kmesh_summary.json").write_text(json.dumps(dict(recommended=rec, rows=rows), indent=1))


if __name__ == "__main__":
    main()
