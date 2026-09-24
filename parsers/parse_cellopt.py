#!/usr/bin/env python3
"""
parse_cellopt.py -- step 03: converged bulk -> <project>_relaxed.cif + bulk.json.

bulk.json carries what the surface-energy parser needs (E per formula unit
at the production settings) and what 04/05/07 need (relaxed structure,
k-density). Refuses unconverged runs.

  python ~/work_cp2k/parsers/parse_cellopt.py --run-dir 03_cellopt
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
    p.add_argument("--run-dir", default=".")
    a = p.parse_args()

    d = Path(a.run_dir)
    out = co.find_output(d)
    t = co.read(out) if out else ""
    if not co.opt_converged(t):
        sys.exit(f"ERROR: {d}: CELL_OPT not converged (no OPTIMIZATION COMPLETED).")
    rst = co.find_restart(d)
    st = cb.read_structure(rst)
    project = rst.name.replace("-1.restart", "")
    cif = d / f"{project}_relaxed.cif"
    st.to(filename=str(cif))
    nfu = cb.n_formula_units(st)
    e = co.final_energy_ev(t)
    run_meta = json.loads((d / "run_meta.json").read_text()) if (d / "run_meta.json").exists() else {}
    bulk = dict(relaxed_cif=str(cif.resolve()), formula=st.composition.formula, natoms=len(st),
                n_formula_units=nfu, E_total_eV=e, E_per_fu_eV=e / nfu, volume_A3=st.volume,
                abc=list(st.lattice.abc), angles=list(st.lattice.angles),
                kpoints=run_meta.get("kpoints"), k_density=run_meta.get("k_density"),
                cutoff=run_meta.get("cutoff"), rel_cutoff=run_meta.get("rel_cutoff"),
                basis=run_meta.get("basis"), mpi_ranks=co.mpi_ranks(t))
    (d / "bulk.json").write_text(json.dumps(bulk, indent=1))
    print(json.dumps(bulk, indent=1))


if __name__ == "__main__":
    main()
