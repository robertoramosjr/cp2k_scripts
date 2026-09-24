#!/usr/bin/env python3
"""
cutoff_convergence.py

Generates a series of single-point (ENERGY_FORCE) CP2K inputs across a range
of MGRID CUTOFF values (REL_CUTOFF and k-mesh held fixed), for plane-wave-grid
cutoff convergence testing. This should be run BEFORE kpoint_convergence.py:
converge CUTOFF first (cheap, saturates fast), then converge the k-mesh with
CUTOFF/REL_CUTOFF fixed at the converged value -- same dependency order as a
VASP ENCUT scan followed by a KPOINTS scan.

The k-mesh is deliberately kept coarse (Gamma-only by default) during this
scan to keep it cheap -- CUTOFF affects the real-space integration grid and
is largely decoupled from k-point sampling.

Usage:
    python ~/scripts/cutoff_convergence.py --mp-id mp-390 \\
        --project-name TiO2_anatase --cutoff-min 200 --cutoff-max 800 --cutoff-step 100

    python ~/scripts/cutoff_convergence.py --cif TiO2_anatase.cif \\
        --project-name TiO2_anatase --cutoff-list "200,300,400,500,600,800,1000"
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import input_generator as ig  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a MGRID CUTOFF convergence scan of CP2K single-point inputs."
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mp-id", type=str, help="Materials Project ID, e.g. mp-390")
    source.add_argument("--cif", type=str, help="Path to a local structure file (CIF, POSCAR, etc.)")

    parser.add_argument("--api-key", type=str, default=None,
                         help="Materials Project API key. Defaults to $MP_API_KEY env var.")
    parser.add_argument("--project-name", type=str, default="cutoffconv",
                         help="Base PROJECT name; each cutoff appends its own tag.")
    parser.add_argument("--functional", type=str, default="PBE",
                         choices=["PBE", "PBE0", "HSE06"])
    parser.add_argument("--basis-family", type=str, default=ig.BASIS_FAMILY_DEFAULT)
    parser.add_argument("--rel-cutoff", type=int, default=60,
                         help="MGRID REL_CUTOFF, held fixed across the scan.")
    parser.add_argument("--kpoints", type=int, nargs=3, default=[1, 1, 1],
                         metavar=("KX", "KY", "KZ"),
                         help="Fixed Monkhorst-Pack mesh during the cutoff scan "
                              "(deliberately coarse/Gamma by default to keep it cheap).")
    parser.add_argument("--basis-set-file", type=str, default="BASIS_MOLOPT")
    parser.add_argument("--potential-file", type=str, default="GTH_POTENTIALS")
    parser.add_argument("--kind-overrides", type=str, default=None)

    # Convergence knobs -- same VASP-INCAR-derived defaults as input_generator.py
    # (EDIFF, NELM, AMIX). No --max-iter/--max-force here: this scan is
    # ENERGY_FORCE (single-point), there's no geometry/cell optimizer running.
    parser.add_argument("--eps-scf", type=float, default=1.0E-6,
                         help="SCF energy convergence criterion (mirrors VASP's EDIFF).")
    parser.add_argument("--max-scf", type=int, default=300,
                         help="Max inner SCF iterations (loose safety-margin analog of VASP's NELM).")
    parser.add_argument("--outer-max-scf", type=int, default=50,
                         help="Max OUTER_SCF iterations.")
    parser.add_argument("--mixing-alpha", type=float, default=0.4,
                         help="Density mixing fraction for non-Gamma k-meshes (mirrors VASP's AMIX).")

    parser.add_argument("--cutoff-min", type=int, default=200, help="Smallest CUTOFF (Ry).")
    parser.add_argument("--cutoff-max", type=int, default=1000, help="Largest CUTOFF (Ry).")
    parser.add_argument("--cutoff-step", type=int, default=100, help="Step between CUTOFF values (Ry).")
    parser.add_argument("--cutoff-list", type=str, default=None,
                         help='Explicit CUTOFF values, e.g. "200,300,400,600,800". '
                              'Overrides --cutoff-min/--cutoff-max/--cutoff-step.')

    parser.add_argument("--output-dir", type=str, default="cutoff_convergence",
                         help="Root directory where per-cutoff subfolders are created.")

    return parser.parse_args()


def build_cutoff_list(args):
    if args.cutoff_list:
        return [int(x) for x in args.cutoff_list.split(",")]
    return list(range(args.cutoff_min, args.cutoff_max + 1, args.cutoff_step))


def render_scan_input(cell_block, coord_block, kind_blocks, xc_block, args,
                       project_name, cutoff):
    kx, ky, kz = args.kpoints
    return f"""&GLOBAL
  PROJECT  {project_name}
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME  {args.basis_set_file}
    POTENTIAL_FILE_NAME  {args.potential_file}
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
    &END QS
    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF {args.rel_cutoff}
      NGRIDS 5
    &END MGRID
    {xc_block}
    {ig.build_scf_block([kx, ky, kz], disable_wfn_restart=True,
                         eps_scf=args.eps_scf, max_scf=args.max_scf,
                         outer_max_scf=args.outer_max_scf, mixing_alpha=args.mixing_alpha)}
{ig.build_kpoints_block([kx, ky, kz])}  &END DFT
  &SUBSYS
    {cell_block}
    {coord_block}
    {kind_blocks}
  &END SUBSYS
&END FORCE_EVAL
"""


def main():
    args = parse_args()
    structure = ig.get_structure(mp_id=args.mp_id, cif_path=args.cif, api_key=args.api_key)
    overrides = ig.load_kind_overrides(args.kind_overrides)

    cell_block = ig.build_cell_block(structure)
    coord_block = ig.build_coord_block(structure)
    kind_blocks = ig.build_kind_blocks(
        structure, args.basis_family, args.functional, overrides, use_admm=False
    )
    xc_block = ig.build_xc_block(args.functional)

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    cutoffs = build_cutoff_list(args)
    print(f"Structure: {structure.composition.reduced_formula} ({len(structure)} sites)")
    print(f"Fixed k-mesh for this scan: {args.kpoints[0]}x{args.kpoints[1]}x{args.kpoints[2]}, "
          f"REL_CUTOFF={args.rel_cutoff}")
    print(f"Generating {len(cutoffs)} cutoff convergence inputs in '{out_root}/'\n")

    for cutoff in cutoffs:
        tag = f"cutoff{cutoff}"
        run_dir = out_root / tag
        run_dir.mkdir(exist_ok=True)
        project_name = f"{args.project_name}_{tag}"

        input_text = render_scan_input(
            cell_block, coord_block, kind_blocks, xc_block, args, project_name, cutoff
        )
        input_path = run_dir / f"{project_name}.inp"
        input_path.write_text(input_text, encoding="utf-8")
        print(f"  {tag:<14} -> {input_path}")

    print(f"\nAtoms in cell (needed later for meV/atom normalization): {len(structure)}")
    print("Run each with, e.g.:")
    print(f'  for d in {out_root}/*/; do (cd "$d" && cp2k.psmp -i *.inp -o *.inp.out); done')
    print("\nThen parse with:")
    print(f"  python ~/scripts/parse_cutoff_convergence.py --scan-dir {out_root} "
          f"--n-atoms {len(structure)}")
    print("\nOnce CUTOFF is converged, feed it into kpoint_convergence.py via --cutoff "
          "to run the k-mesh scan on top of it.")


if __name__ == "__main__":
    main()