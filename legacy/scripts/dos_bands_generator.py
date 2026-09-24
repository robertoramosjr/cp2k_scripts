#!/usr/bin/env python3
"""
dos_bands_generator.py

Generates a single CP2K input that computes BOTH the density of states
(PDOS) and the electronic band structure along the crystal's standard
high-symmetry k-path, from one converged SCF on the relaxed structure.

Why this works in one run: both DOS and band structure are evaluated from
the SAME self-consistent ground-state density obtained from a single SCF
on the production Monkhorst-Pack mesh (RUN_TYPE ENERGY, no relaxation).
CP2K's &DFT %PRINT%PDOS uses that density directly; &DFT%PRINT%BAND_STRUCTURE
re-diagonalizes the converged Hamiltonian along an explicit k-point path
(non-self-consistently) -- no second SCF needed.

The high-symmetry path is generated automatically via pymatgen's
HighSymmKpath (Setyawan-Curtarolo convention), which also returns the
STANDARDIZED PRIMITIVE cell the path's fractional coordinates are defined
against. That standardized cell -- not your original relaxed structure's
cell/orientation -- is what gets written into &CELL/&COORD, because the
k-path coordinates are only meaningful relative to it.

Usage:
    python ~/scripts/dos_bands_generator.py --cif TiO2_etot_relaxed.cif \\
        --project-name TiO2_dos_bands --cutoff 600 --kpoints 6 6 6 \\
        --added-mos 20 --wfn-restart TiO2_etot-RESTART.wfn
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import input_generator as ig  # noqa: E402

try:
    from pymatgen.symmetry.bandstructure import HighSymmKpath
except ImportError:
    sys.exit("ERROR: pymatgen is required. Install with: python3 -m pip install pymatgen")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a combined DOS + band-structure CP2K input."
    )
    parser.add_argument("--cif", type=str, required=True,
                         help="Path to the RELAXED structure (e.g. from extract_relaxed_structure.py).")
    parser.add_argument("--project-name", type=str, default="dos_bands")
    parser.add_argument("--functional", type=str, default="PBE",
                         choices=["PBE", "PBE0", "HSE06"])
    parser.add_argument("--basis-family", type=str, default=ig.BASIS_FAMILY_DEFAULT)
    parser.add_argument("--cutoff", type=int, required=True,
                         help="MGRID CUTOFF (Ry), from the converged scan.")
    parser.add_argument("--rel-cutoff", type=int, default=60)
    parser.add_argument("--kpoints", type=int, nargs=3, required=True,
                         metavar=("KX", "KY", "KZ"),
                         help="Etot's production Monkhorst-Pack mesh (baseline, pre-densification).")
    parser.add_argument("--dos-kpoints", type=int, nargs=3, default=None,
                         metavar=("KX", "KY", "KZ"),
                         help="Explicit denser mesh for this DOS/bands SCF. Overrides "
                              "--dos-kpoints-multiplier if given.")
    parser.add_argument("--dos-kpoints-multiplier", type=float, default=2.0,
                         help="If --dos-kpoints is not given, multiply --kpoints by this "
                              "factor (rounded) to get the DOS/bands mesh. DOS/PDOS need a "
                              "denser mesh than relaxation-stage energies/forces do -- same "
                              "reason you multiply KPOINTS by 2-3x for DOS in VASP. Note: "
                              "the band-structure PATH itself doesn't need this (CP2K "
                              "Fourier-interpolates it from real-space KS matrices), but the "
                              "PDOS is tied directly to this SCF mesh, so it does.")
    parser.add_argument("--basis-set-file", type=str, default="BASIS_MOLOPT")
    parser.add_argument("--potential-file", type=str, default="GTH_POTENTIALS")
    parser.add_argument("--kind-overrides", type=str, default=None)

    parser.add_argument("--eps-scf", type=float, default=1.0E-6)
    parser.add_argument("--max-scf", type=int, default=300)
    parser.add_argument("--outer-max-scf", type=int, default=50)
    parser.add_argument("--added-mos", type=int, default=None,
                         help="Extra empty bands in the main SCF, so PDOS shows the "
                              "conduction band. Defaults to auto: total occupied MOs "
                              "(valence electrons / 2), doubling total bands -- the CP2K "
                              "analog of VASP's NBANDS=NELECT.")
    parser.add_argument("--wfn-restart", type=str, default=None,
                         help="Optional path to a converged .wfn (e.g. from Etot) to "
                              "seed SCF_GUESS RESTART instead of ATOMIC.")

    parser.add_argument("--band-points", type=int, default=20,
                         help="Number of k-points sampled along each path segment.")
    parser.add_argument("--band-added-mos", type=int, default=None,
                         help="ADDED_MOS specific to the BAND_STRUCTURE print block. "
                              "Defaults to the same value as --added-mos if not given.")

    parser.add_argument("--verify-relax", action="store_true",
                         help="Instead of writing the DOS/bands input, write a quick "
                              "confirmation GEO_OPT input on the pymatgen-standardized "
                              "primitive cell (used for the band path). Run this FIRST: "
                              "pymatgen's symmetrization can nudge atomic positions "
                              "slightly relative to your true relaxed minimum, so this "
                              "confirms (cheaply -- should converge in ~1 step if it's "
                              "really just a rigid relabeling) that the standardized cell "
                              "is still at the energy minimum before computing DOS/bands "
                              "on it. Feed the result back in as --cif for the real run.")

    parser.add_argument("--output", type=str, default=None,
                         help="Output .inp path. Defaults to <project-name>.inp.")
    return parser.parse_args()


def sanitize_label(label: str) -> str:
    """CP2K SPECIAL_POINT names can't contain LaTeX-style backslashes etc."""
    return label.replace("\\", "").replace("$", "").replace(" ", "").upper() or "PT"


def build_kpoint_set_blocks(kpath, band_points: int) -> str:
    branches = kpath.kpath["path"]
    points = kpath.kpath["kpoints"]

    blocks = []
    for branch in branches:
        for start_label, end_label in zip(branch[:-1], branch[1:]):
            start_coord = points[start_label]
            end_coord = points[end_label]
            start_name = sanitize_label(start_label)
            end_name = sanitize_label(end_label)
            blocks.append(
                "      &KPOINT_SET\n"
                f"        NPOINTS {band_points}\n"
                f"        SPECIAL_POINT {start_name}  {start_coord[0]:.8f} {start_coord[1]:.8f} {start_coord[2]:.8f}\n"
                f"        SPECIAL_POINT {end_name}  {end_coord[0]:.8f} {end_coord[1]:.8f} {end_coord[2]:.8f}\n"
                "      &END KPOINT_SET"
            )
    return "\n".join(blocks)


def main():
    args = parse_args()

    cif_path = Path(args.cif)
    if not cif_path.exists():
        sys.exit(f"ERROR: structure file not found: {cif_path}")

    from pymatgen.core import Structure
    relaxed_structure = Structure.from_file(str(cif_path))

    # Standard high-symmetry path (Setyawan-Curtarolo). kpath.prim is the
    # STANDARDIZED PRIMITIVE cell the path's fractional coordinates are
    # defined against -- we MUST simulate this cell, not the original one,
    # or the k-path coordinates won't correspond to the actual lattice.
    kpath = HighSymmKpath(relaxed_structure)
    structure = kpath.prim

    if len(structure) != len(relaxed_structure):
        print(f"NOTE: pymatgen standardized the cell for the k-path "
              f"({len(relaxed_structure)} -> {len(structure)} sites). "
              f"This is expected (primitive-cell reduction) and the input "
              f"below uses the standardized cell consistently.")

    overrides = ig.load_kind_overrides(args.kind_overrides)

    # --- verify-relax mode: just confirm the standardized cell is still at
    # the energy minimum (cheap GEO_OPT), don't build the DOS/bands input yet ---
    if args.verify_relax:
        verify_project = f"{args.project_name}_verify_relax"
        verify_args = argparse.Namespace(
            project_name=verify_project,
            run_type="GEO_OPT",
            functional=args.functional,
            basis_family=args.basis_family,
            cutoff=args.cutoff,
            rel_cutoff=args.rel_cutoff,
            stress_rel_cutoff_factor=1.0,  # irrelevant for GEO_OPT, no stress margin applied
            eps_scf=args.eps_scf,
            max_scf=args.max_scf,
            outer_max_scf=args.outer_max_scf,
            max_iter=20,   # should converge almost immediately if it's just a relabeling
            max_force=0.010,
            kpoints=args.kpoints,
            basis_set_file=args.basis_set_file,
            potential_file=args.potential_file,
            kind_overrides=args.kind_overrides,
        )
        input_text = ig.render_input(structure, verify_args)
        output_path = Path(args.output) if args.output else Path.cwd() / f"{verify_project}.inp"
        output_path.write_text(input_text, encoding="utf-8")
        print(f"Structure (standardized primitive): {structure.composition.reduced_formula} "
              f"({len(structure)} sites)")
        print(f"VERIFY-RELAX input written to: {output_path}")
        print("Run this first. It should converge in ~1 step if the standardization was "
              "just a rigid relabeling of your already-relaxed cell. Extract the result with "
              "extract_relaxed_structure.py and feed THAT back in as --cif for the real DOS/bands run "
              "(without --verify-relax).")
        return

    # --- Auto-size ADDED_MOS from total valence electrons (like VASP's
    # NBANDS=NELECT) unless the user gave an explicit value ---
    if args.added_mos is None:
        n_valence = ig.total_valence_electrons(structure, overrides)
        occupied_mos = -(-n_valence // 2)  # ceil division, non-spin-polarized
        added_mos = occupied_mos
        print(f"Auto ADDED_MOS: {n_valence} valence electrons -> {occupied_mos} occupied MOs "
              f"-> ADDED_MOS={added_mos} (doubles total bands, like VASP's NBANDS=NELECT)")
    else:
        added_mos = args.added_mos

    # --- Denser k-mesh for the DOS/bands SCF (PDOS is tied to this mesh;
    # the band PATH itself doesn't need it, see --dos-kpoints-multiplier help) ---
    if args.dos_kpoints:
        dos_kpoints = args.dos_kpoints
    else:
        dos_kpoints = [max(1, round(k * args.dos_kpoints_multiplier)) for k in args.kpoints]
    print(f"DOS/bands SCF k-mesh: {dos_kpoints[0]}x{dos_kpoints[1]}x{dos_kpoints[2]} "
          f"(Etot mesh was {args.kpoints[0]}x{args.kpoints[1]}x{args.kpoints[2]})")

    cell_block = ig.build_cell_block(structure)
    coord_block = ig.build_coord_block(structure)
    kind_blocks = ig.build_kind_blocks(structure, args.basis_family, args.functional,
                                        overrides, use_admm=(args.functional in ("PBE0", "HSE06")))
    xc_block = ig.build_xc_block(args.functional)
    scf_block = ig.build_scf_block(
        dos_kpoints, disable_wfn_restart=True,  # one-shot ENERGY run, no continuation needed
        eps_scf=args.eps_scf, max_scf=args.max_scf, outer_max_scf=args.outer_max_scf,
        added_mos=added_mos, wfn_restart_file=args.wfn_restart,
    )

    kpoint_set_blocks = build_kpoint_set_blocks(kpath, args.band_points)
    band_added_mos = args.band_added_mos if args.band_added_mos is not None else added_mos
    kx, ky, kz = dos_kpoints

    template = f"""&GLOBAL
  PROJECT  {args.project_name}
  RUN_TYPE ENERGY
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
      CUTOFF {args.cutoff}
      REL_CUTOFF {args.rel_cutoff}
      NGRIDS 5
    &END MGRID
    {xc_block}
    {scf_block}
    &KPOINTS
      SCHEME MONKHORST-PACK {kx} {ky} {kz}
    &END KPOINTS
    &PRINT
      &PDOS
        COMPONENTS
        NLUMO -1
      &END PDOS
      &BAND_STRUCTURE
        ADDED_MOS {band_added_mos}
        FILE_NAME {args.project_name}.bs
{kpoint_set_blocks}
      &END BAND_STRUCTURE
    &END PRINT
  &END DFT
  &SUBSYS
    {cell_block}
    {coord_block}
    {kind_blocks}
  &END SUBSYS
&END FORCE_EVAL
"""

    output_path = Path(args.output) if args.output else Path.cwd() / f"{args.project_name}.inp"
    output_path.write_text(template, encoding="utf-8")

    print(f"Structure (standardized primitive): {structure.composition.reduced_formula} "
          f"({len(structure)} sites)")
    print(f"Band path: {' -> '.join(' -> '.join(sanitize_label(l) for l in b) for b in kpath.kpath['path'])}")
    print(f"CP2K input written to: {output_path}")
    print(f"PDOS files: <project>-k*-*.pdos ; band structure: {args.project_name}.bs")


if __name__ == "__main__":
    main()