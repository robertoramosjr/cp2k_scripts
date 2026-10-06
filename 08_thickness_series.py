#!/usr/bin/env python3
"""
08_thickness_series.py -- slab-thickness convergence series for ONE production slab.

Re-cuts the plane of --slab-dir (a folder from 06) with 06_slab_cut.py at each
--min-slab value, picks in every cut the SAME termination as the production
slab (matched by broken bonds per area and construction, not by index: the
ranking index can change with thickness), and writes a 07 input for it in
<output-dir>/<slab>_thk<T>/ with the production settings (bulk rescaling,
k-density, vacuum, CUTOFF/REL_CUTOFF ... passed through to 07 unchanged).
Cuts that give the same number of atoms as a previous one are skipped.

  python ~/work_cp2k/08_thickness_series.py --slab-dir slabs_studies/slab_-111_1 \\
      --min-slabs 8 12 16 20 --output-dir thickness_convergence \\
      --bulk-cp2k 03_cellopt/ZrO2_m_k3-3-3_relaxed.cif --k-density 16.12 \\
      --no-dipole --vacuum 25 --basis DZVP-MOLOPT-SR-GTH --cutoff 600 --rel-cutoff 50
Then run the folders (job_sequential.sh / job_scan_array.sh) and
  python ~/work_cp2k/parsers/parse_thickness.py --scan-dir thickness_convergence \\
      --bulk-json 03_cellopt/bulk.json
Default run type is GEO_OPT (gamma of relaxed slabs is what converges with
thickness); pass --run-type ENERGY for a cheaper unrelaxed check.
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slab-dir", required=True, help="Production slab folder from 06 (slab_meta.json).")
    p.add_argument("--min-slabs", type=float, nargs="+", required=True, help="06 --min-slab values (angstrom).")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--tol-broken", type=float, default=1e-3,
                   help="Match tolerance on broken bonds per A^2 (default 1e-3).")
    a, passthrough = p.parse_known_args()   # everything else goes to 07_slab_opt.py

    src = Path(a.slab_dir).resolve()
    meta = json.loads((src / "slab_meta.json").read_text())
    bulk = Path(meta["bulk"])
    if not bulk.exists():  # path recorded on another machine: use this phase's copy
        bulk = src.parents[1] / bulk.name
    hkl = ",".join(str(h) for h in meta["hkl"])
    out_root = Path(a.output_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    seen_natoms = set()

    for t in sorted(a.min_slabs):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run([sys.executable, str(HERE / "06_slab_cut.py"), "--bulk", str(bulk),
                                f"--planes={hkl}", "--min-slab", str(t), "--vacuum", str(meta["vacuum_A"]),
                                "--max-terminations", "20", "--bond-cutoff", str(meta["bond_cutoff_A"]),
                                "--output-dir", tmp], capture_output=True, text=True)
            if r.returncode:
                sys.exit(f"06_slab_cut.py failed at --min-slab {t}:\n{r.stdout}\n{r.stderr}")
            match = None
            for d in sorted(Path(tmp).glob("slab_*")):
                m = json.loads((d / "slab_meta.json").read_text())
                if (m["origin"] == meta["origin"]
                        and abs(m["broken_bonds_per_A2"] - meta["broken_bonds_per_A2"]) < a.tol_broken
                        and abs(m["area_A2"] - meta["area_A2"]) < 1e-3 * meta["area_A2"]):
                    match = (d, m)
                    break
            if match is None:
                print(f"[WARNING] min-slab {t:g}: termination of {src.name} not found in the cut, skipped")
                continue
            d, m = match
            if m["natoms"] in seen_natoms:
                print(f"[INFO] min-slab {t:g}: {m['natoms']} at, same slab as a thinner value, skipped")
                continue
            seen_natoms.add(m["natoms"])
            name = f"{src.name}_thk{t:g}"
            dest = out_root / name
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(d, dest)
            m["bulk"] = str(bulk)
            m["termination"] = meta["termination"]          # production label, not this cut's index
            m["thickness_series_of"] = str(src)
            (dest / "slab_meta.json").write_text(json.dumps(m, indent=1))
        r = subprocess.run([sys.executable, str(HERE / "07_slab_opt.py"), "--slab-dir", str(dest),
                            "--project", name, *passthrough], capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"07_slab_opt.py failed for {name}:\n{r.stdout}\n{r.stderr}")
        print(f"min-slab {t:5g}: {m['natoms']:4d} at, t={m['thickness_A']:6.2f} A -> {dest}\n   {r.stdout.strip()}")


if __name__ == "__main__":
    main()
