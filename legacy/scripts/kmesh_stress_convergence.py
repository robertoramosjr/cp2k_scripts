#!/usr/bin/env python3
"""
kmesh_stress_convergence.py

Generates a series of FULL CELL_OPT ("stress") CP2K inputs across a range of
Monkhorst-Pack k-point meshes. This is the rigorous k-convergence methodology
(mirrors the user's VASP "stress" job): instead of testing raw SCF energy on
a fixed unrelaxed geometry, each k-mesh point is relaxed to its own
converged cell -- and the RELAXED energy/volume are what get compared across
meshes. More expensive than kpoint_convergence.py, but directly tests the
quantity you actually care about downstream.

Order of operations (same dependency chain as ENCUT->KPOINTS in VASP):
    1. cutoff_convergence.py   -- cheap, single-point CUTOFF scan
    2. kmesh_stress_convergence.py (this script) -- full CELL_OPT at each
       k-mesh, using the converged CUTOFF and a REL_CUTOFF stress-safety
       margin (see input_generator.py --stress-rel-cutoff-factor)
    3. parse_kmesh_stress_convergence.py -- pick the converged k-mesh from
       relaxed energy AND relaxed volume
    4. A single production GEO_OPT ("Etot") at that k-mesh, standard
       REL_CUTOFF, starting from the relaxed structure -- see
       input_generator.py --run-type GEO_OPT

Usage:
    python ~/scripts/kmesh_stress_convergence.py --mp-id mp-390 \\
        --project-name TiO2_anatase --cutoff 600 --k-min 2 --k-max 8 --k-step 2
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import input_generator as ig  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a k-mesh convergence scan of FULL CP2K CELL_OPT (stress) inputs."
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mp-id", type=str, help="Materials Project ID, e.g. mp-390")
    source.add_argument("--cif", type=str, help="Path to a local structure file (CIF, POSCAR, etc.)")

    parser.add_argument("--api-key", type=str, default=None,
                         help="Materials Project API key. Defaults to $MP_API_KEY env var.")
    parser.add_argument("--project-name", type=str, default="stress",
                         help="Base PROJECT name; each mesh appends its own tag.")
    parser.add_argument("--functional", type=str, default="PBE",
                         choices=["PBE", "PBE0", "HSE06"])
    parser.add_argument("--basis-family", type=str, default=ig.BASIS_FAMILY_DEFAULT)
    parser.add_argument("--cutoff", type=int, required=True,
                         help="MGRID CUTOFF (Ry), already converged via cutoff_convergence.py.")
    parser.add_argument("--rel-cutoff", type=int, default=60,
                         help="Base MGRID REL_CUTOFF (Ry) before the stress safety margin.")
    parser.add_argument("--stress-rel-cutoff-factor", type=float, default=1.5,
                         help="Safety multiplier on REL_CUTOFF for the stress tensor "
                              "(applied automatically by input_generator.render_input "
                              "since RUN_TYPE=CELL_OPT). Set to 1.0 to disable.")
    parser.add_argument("--basis-set-file", type=str, default="BASIS_MOLOPT")
    parser.add_argument("--potential-file", type=str, default="GTH_POTENTIALS")
    parser.add_argument("--kind-overrides", type=str, default=None)

    # Convergence / optimizer knobs -- same VASP-INCAR-derived defaults as
    # input_generator.py (EDIFF, NELM, NSW, EDIFFG). See build_scf_block's
    # docstring for why NELM isn't copied 1:1 into max_scf.
    parser.add_argument("--eps-scf", type=float, default=1.0E-6)
    parser.add_argument("--max-scf", type=int, default=300)
    parser.add_argument("--outer-max-scf", type=int, default=50)
    parser.add_argument("--mixing-alpha", type=float, default=0.4,
                         help="Density mixing fraction for non-Gamma k-meshes (mirrors VASP's AMIX).")
    parser.add_argument("--max-iter", type=int, default=500,
                         help="Max CELL_OPT optimizer steps (mirrors VASP's NSW).")
    parser.add_argument("--max-force", type=float, default=0.010,
                         help="Force convergence in eV/angstrom (mirrors VASP's EDIFFG).")

    parser.add_argument("--k-min", type=int, default=2, help="Smallest isotropic mesh N (NxNxN).")
    parser.add_argument("--k-max", type=int, default=8, help="Largest isotropic mesh N (NxNxN).")
    parser.add_argument("--k-step", type=int, default=2, help="Step between isotropic meshes.")
    parser.add_argument("--k-list", type=str, default=None,
                         help='Explicit meshes, e.g. "2,2,2;4,4,4;6,6,6". '
                              'Overrides --k-min/--k-max/--k-step.')

    parser.add_argument("--output-dir", type=str, default="kmesh_stress_convergence",
                         help="Root directory where per-mesh subfolders are created.")

    return parser.parse_args()


def build_mesh_list(args):
    if args.k_list:
        meshes = []
        for triple in args.k_list.split(";"):
            kx, ky, kz = (int(x) for x in triple.split(","))
            meshes.append((kx, ky, kz))
        return meshes
    return [(n, n, n) for n in range(args.k_min, args.k_max + 1, args.k_step)]


def main():
    args = parse_args()
    structure = ig.get_structure(mp_id=args.mp_id, cif_path=args.cif, api_key=args.api_key)

    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    meshes = build_mesh_list(args)
    print(f"Structure: {structure.composition.reduced_formula} ({len(structure)} sites)")
    print(f"Generating {len(meshes)} full CELL_OPT (stress) inputs in '{out_root}/'")
    print("(each one relaxes cell + ions at its own k-mesh -- this is the expensive, "
          "rigorous methodology, not a single-point scan)\n")

    for kx, ky, kz in meshes:
        tag = f"k{kx}-{ky}-{kz}"
        run_dir = out_root / tag
        run_dir.mkdir(exist_ok=True)
        project_name = f"{args.project_name}_{tag}"

        # Reuse input_generator's render_input directly -- avoids duplicating the
        # FORCE_EVAL/DFT/SUBSYS template (and its STRESS_TENSOR / REL_CUTOFF logic)
        # in a second place.
        run_args = argparse.Namespace(
            project_name=project_name,
            run_type="CELL_OPT",
            functional=args.functional,
            basis_family=args.basis_family,
            cutoff=args.cutoff,
            rel_cutoff=args.rel_cutoff,
            stress_rel_cutoff_factor=args.stress_rel_cutoff_factor,
            eps_scf=args.eps_scf,
            max_scf=args.max_scf,
            outer_max_scf=args.outer_max_scf,
            mixing_alpha=args.mixing_alpha,
            max_iter=args.max_iter,
            max_force=args.max_force,
            kpoints=[kx, ky, kz],
            basis_set_file=args.basis_set_file,
            potential_file=args.potential_file,
            kind_overrides=args.kind_overrides,
        )
        input_text = ig.render_input(structure, run_args)
        input_path = run_dir / f"{project_name}.inp"
        input_path.write_text(input_text, encoding="utf-8")
        print(f"  {tag:<12} -> {input_path}")

    print(f"\nAtoms in cell (needed later for meV/atom normalization): {len(structure)}")
    print("Run each with:")
    print(f"  sbatch ~/scripts/job_kmesh_stress_convergence.sh {out_root}")
    print("\nThen parse with:")
    print(f"  python ~/scripts/parse_kmesh_stress_convergence.py --scan-dir {out_root} "
          f"--n-atoms {len(structure)}")


if __name__ == "__main__":
    main()