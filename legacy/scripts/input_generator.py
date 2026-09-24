#!/usr/bin/env python3
"""
input_generator.py

Generates a CP2K input file (.inp) from a Materials Project structure
(by material_id) or a local structure file (CIF, POSCAR, etc.), using
pymatgen for structure parsing.

Intended usage pattern (run from the working directory of the project):

    python ~/scripts/input_generator.py --mp-id mp-390 --run-type CELL_OPT
    python ~/scripts/input_generator.py --cif TiO2_anatase.cif --functional HSE06

Requirements:
    pip install pymatgen mp-api

Materials Project API key:
    Set the MP_API_KEY environment variable, or pass --api-key explicitly.
    Get a key at https://next-gen.materialsproject.org/api

Author: generated for Roberto de Aguiar Ramos Junior
"""

import argparse
import json
import os
import sys
from pathlib import Path

try:
    from pymatgen.core import Structure
except ImportError:
    sys.exit("ERROR: pymatgen is required. Install with: pip install pymatgen")


# ----------------------------------------------------------------------
# Default GTH pseudopotential valence-electron counts (the "qN" suffix)
# for elements commonly used in oxide / perovskite DFT work.
# These follow the standard GTH_POTENTIALS / BASIS_MOLOPT naming
# convention shipped with CP2K. ALWAYS cross-check against your local
# GTH_POTENTIALS and BASIS_MOLOPT files before production runs, since
# library versions can differ.
# ----------------------------------------------------------------------
DEFAULT_GTH_VALENCE = {
    "H": 1, "Li": 3, "Be": 4, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7,
    "Na": 9, "Mg": 10, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7,
    "K": 9, "Ca": 10, "Sc": 11, "Ti": 12, "V": 13, "Cr": 14, "Mn": 15,
    "Fe": 16, "Co": 17, "Ni": 18, "Cu": 11, "Zn": 12,
    "Sr": 10, "Y": 11, "Zr": 12, "Nb": 13, "Mo": 14,
    "Sn": 4, "Ba": 10, "Hf": 12, "Pb": 4,
}

BASIS_FAMILY_DEFAULT = "MOLOPT-SR"  # short-range MOLOPT basis, robust for periodic solids

def build_motion_block(run_type: str, max_iter: int, max_force_ev_ang: float) -> str:
    """Build the &MOTION block, with MAX_ITER/MAX_FORCE exposed as CLI-tunable
    knobs (mirrors VASP's NSW and EDIFFG). MAX_FORCE is written directly in
    eV/angstrom using CP2K's compound-unit syntax, matching EDIFFG's units
    with no manual Hartree/Bohr conversion needed.
    """
    if run_type == "CELL_OPT":
        return (
            "&MOTION\n"
            "  &CELL_OPT\n"
            "    OPTIMIZER BFGS\n"
            f"    MAX_ITER {max_iter}\n"
            f"    MAX_FORCE [eV*angstrom^-1] {max_force_ev_ang}\n"
            "    PRESSURE_TOLERANCE [bar] 100\n"
            "  &END CELL_OPT\n"
            "  &PRINT\n"
            "    &RESTART\n"
            "      BACKUP_COPIES 3\n"
            "      &EACH\n"
            "        CELL_OPT 1\n"
            "      &END EACH\n"
            "    &END RESTART\n"
            "  &END PRINT\n"
            "&END MOTION\n"
        )
    if run_type == "GEO_OPT":
        return (
            "&MOTION\n"
            "  &GEO_OPT\n"
            "    OPTIMIZER BFGS\n"
            f"    MAX_ITER {max_iter}\n"
            f"    MAX_FORCE [eV*angstrom^-1] {max_force_ev_ang}\n"
            "  &END GEO_OPT\n"
            "  &PRINT\n"
            "    &RESTART\n"
            "      BACKUP_COPIES 3\n"
            "      &EACH\n"
            "        GEO_OPT 1\n"
            "      &END EACH\n"
            "    &END RESTART\n"
            "  &END PRINT\n"
            "&END MOTION\n"
        )
    return ""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a CP2K input file from a Materials Project ID or a local structure file."
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mp-id", type=str, help="Materials Project ID, e.g. mp-390")
    source.add_argument("--cif", type=str, help="Path to a local structure file (CIF, POSCAR, etc.)")

    parser.add_argument("--api-key", type=str, default=None,
                         help="Materials Project API key. Defaults to $MP_API_KEY env var.")
    parser.add_argument("--project-name", type=str, default="cp2k_project",
                         help="PROJECT name in &GLOBAL (used as output filename stem).")
    parser.add_argument("--run-type", type=str, default="CELL_OPT",
                         choices=["CELL_OPT", "GEO_OPT", "ENERGY", "ENERGY_FORCE"],
                         help="RUN_TYPE / MOTION block to generate.")
    parser.add_argument("--functional", type=str, default="PBE",
                         choices=["PBE", "PBE0", "HSE06"],
                         help="Exchange-correlation functional.")
    parser.add_argument("--basis-family", type=str, default=BASIS_FAMILY_DEFAULT,
                         help="MOLOPT basis family, e.g. MOLOPT-SR or MOLOPT.")
    parser.add_argument("--cutoff", type=int, default=500, help="MGRID CUTOFF in Ry.")
    parser.add_argument("--rel-cutoff", type=int, default=60, help="MGRID REL_CUTOFF in Ry.")
    parser.add_argument("--stress-rel-cutoff-factor", type=float, default=1.5,
                         help="Safety multiplier applied to REL_CUTOFF ONLY for RUN_TYPE=CELL_OPT. "
                              "CP2K's stress tensor is more sensitive to REL_CUTOFF completeness "
                              "than plain energies/forces are (see manual: MGRID, "
                              "Geometry and cell optimization). This is the GPW analog of "
                              "the VASP practice of using a higher ENCUT for volume relaxation "
                              "-- but here it targets REL_CUTOFF specifically, not CUTOFF, "
                              "since the Gaussian basis itself does not depend on cell volume. "
                              "Set to 1.0 to disable.")
    parser.add_argument("--kpoints", type=int, nargs=3, default=[6, 6, 6],
                         metavar=("KX", "KY", "KZ"), help="Monkhorst-Pack k-point mesh.")
    parser.add_argument("--basis-set-file", type=str, default="BASIS_MOLOPT",
                         help="BASIS_SET_FILE_NAME value.")
    parser.add_argument("--potential-file", type=str, default="GTH_POTENTIALS",
                         help="POTENTIAL_FILE_NAME value.")
    parser.add_argument("--kind-overrides", type=str, default=None,
                         help="Path to a JSON file overriding basis/potential per element, "
                              "e.g. {\"Ti\": {\"basis_set\": \"...\", \"potential\": \"...\"}}")

    # Convergence / optimizer knobs -- defaults mirror a VASP production INCAR
    # (EDIFF=1.0E-6, NELM=1500, NSW=800, EDIFFG=-0.01) where a sensible CP2K
    # analog exists. NELM is NOT copied 1:1 into MAX_SCF -- see build_scf_block.
    parser.add_argument("--eps-scf", type=float, default=1.0E-6,
                         help="SCF energy convergence criterion (mirrors VASP's EDIFF).")
    parser.add_argument("--max-scf", type=int, default=300,
                         help="Max inner SCF iterations (loosely mirrors VASP's NELM, "
                              "but CP2K's OT/DIIS converges much faster than VASP's "
                              "blocked Davidson -- this is a safety margin, not a 1:1 copy).")
    parser.add_argument("--outer-max-scf", type=int, default=50,
                         help="Max OUTER_SCF iterations.")
    parser.add_argument("--mixing-alpha", type=float, default=0.4,
                         help="Density mixing fraction for non-Gamma k-meshes (mirrors VASP's AMIX).")
    parser.add_argument("--surface-dipole", action="store_true",
                         help="Slab runs: add SURFACE_DIPOLE_CORRECTION along Z (mirrors VASP's "
                              "LDIPOL=.TRUE./IDIPOL=3). The slab normal MUST be along z and the "
                              "vacuum along c -- slab_generator.py reorients cells to guarantee it.")
    parser.add_argument("--max-iter", type=int, default=500,
                         help="Max GEO_OPT/CELL_OPT optimizer steps (mirrors VASP's NSW).")
    parser.add_argument("--max-force", type=float, default=0.010,
                         help="Force convergence criterion in eV/angstrom "
                              "(mirrors VASP's EDIFFG, e.g. -0.01 -> 0.010).")

    parser.add_argument("--output", type=str, default=None,
                         help="Output .inp path. Defaults to <project-name>.inp in the current directory.")

    return parser.parse_args()


def get_structure(mp_id: str = None, cif_path: str = None, api_key: str = None) -> Structure:
    """Fetch a pymatgen Structure either from Materials Project (mp_id) or a local file (cif_path).

    Decoupled from argparse on purpose, so other scripts (e.g. kpoint_convergence.py)
    can import and reuse this function directly.
    """
    if mp_id:
        try:
            from mp_api.client import MPRester
        except Exception as e:
            sys.exit(
                f"ERROR: failed to import mp_api.client.MPRester "
                f"({type(e).__name__}: {e})\n"
                f"If 'pip show mp-api' confirms it IS installed in this environment, "
                f"this is most likely a dependency/version conflict (mp-api vs. pymatgen "
                f"vs. emmet-core), not a missing package. Run this directly to see the "
                f"full traceback:\n"
                f'  python3 -c "from mp_api.client import MPRester"'
            )

        api_key = api_key or os.environ.get("MP_API_KEY")
        if not api_key:
            sys.exit("ERROR: no Materials Project API key found. "
                      "Set MP_API_KEY or pass --api-key.")

        with MPRester(api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mp_id)
        return structure

    path = Path(cif_path)
    if not path.exists():
        sys.exit(f"ERROR: structure file not found: {path}")
    return Structure.from_file(str(path))


def load_kind_overrides(path):
    if path is None:
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_cell_block(structure: Structure) -> str:
    matrix = structure.lattice.matrix
    lines = ["&CELL"]
    for label, vec in zip(("A", "B", "C"), matrix):
        lines.append(f"  {label}   {vec[0]:.8f}   {vec[1]:.8f}   {vec[2]:.8f}")
    lines.append("  PERIODIC XYZ")
    lines.append("&END CELL")
    return "\n".join(lines)


def build_coord_block(structure: Structure) -> str:
    lines = ["&COORD", "  SCALED"]
    for site in structure.sites:
        symbol = site.specie.symbol
        fx, fy, fz = site.frac_coords
        lines.append(f"  {symbol:<4}{fx: .8f}  {fy: .8f}  {fz: .8f}")
    lines.append("&END COORD")
    return "\n".join(lines)


def build_kind_blocks(structure: Structure, basis_family: str, functional: str,
                       overrides: dict, use_admm: bool) -> str:
    functional_tag = "PBE" if functional in ("PBE", "PBE0", "HSE06") else functional
    elements = sorted({site.specie.symbol for site in structure.sites})

    blocks = []
    for element in elements:
        override = overrides.get(element, {})

        if "basis_set" in override and "potential" in override:
            basis_set = override["basis_set"]
            potential = override["potential"]
        else:
            valence = DEFAULT_GTH_VALENCE.get(element)
            if valence is None:
                sys.exit(
                    f"ERROR: no default GTH valence known for element '{element}'. "
                    f"Add it to DEFAULT_GTH_VALENCE or supply --kind-overrides."
                )
            basis_set = f"DZVP-{basis_family}-GTH-q{valence}"
            potential = f"GTH-{functional_tag}-q{valence}"

        block = [f"&KIND {element}", f"  BASIS_SET {basis_set}", f"  POTENTIAL {potential}"]
        if use_admm:
            aux_basis = override.get("aux_basis_set", f"cFIT3-{functional_tag}-q{DEFAULT_GTH_VALENCE.get(element, '')}")
            block.append(f"  BASIS_SET AUX_FIT {aux_basis}")
        block.append(f"&END KIND")
        blocks.append("\n".join(block))

    return "\n".join(blocks)


def build_xc_block(functional: str) -> str:
    """
    NOTE: HSE06 / PBE0 blocks below are a reasonable starting point (LIBXC hybrid
    + truncated HF exchange + ADMM), but hybrid-functional band structure runs
    need extra tuning (screening parameter, EPS_SCHWARZ, ADMM purification method)
    that we will refine together in the band-structure step of the workflow.
    """
    if functional == "PBE":
        return (
            "&XC\n"
            "  &XC_FUNCTIONAL PBE\n"
            "  &END XC_FUNCTIONAL\n"
            "&END XC"
        )

    if functional == "PBE0":
        return (
            "&XC\n"
            "  &XC_FUNCTIONAL\n"
            "    &PBE\n"
            "      SCALE_X 0.75\n"
            "    &END PBE\n"
            "    &PBE_HOLE_T_C_LR\n"
            "    &END PBE_HOLE_T_C_LR\n"
            "  &END XC_FUNCTIONAL\n"
            "  &HF\n"
            "    FRACTION 0.25\n"
            "    &INTERACTION_POTENTIAL\n"
            "      POTENTIAL_TYPE TRUNCATED\n"
            "      CUTOFF_RADIUS 6.0\n"
            "    &END INTERACTION_POTENTIAL\n"
            "  &END HF\n"
            "&END XC"
        )

    # HSE06 -- this CP2K build has no generic &LIBXC wrapper section; each
    # libxc functional is its own named subsection under &XC_FUNCTIONAL
    # instead (confirmed via `cp2k.psmp --xml`: &LIBXC doesn't exist, but
    # &HYB_GGA_XC_HSE06 does, with SCALE/_BETA/_OMEGA_HF/_OMEGA_PBE all
    # optional -- defaults already match the literature HSE06 definition).
    return (
        "&XC\n"
        "  &XC_FUNCTIONAL\n"
        "    &HYB_GGA_XC_HSE06\n"
        "    &END HYB_GGA_XC_HSE06\n"
        "  &END XC_FUNCTIONAL\n"
        "  &HF\n"
        "    FRACTION 0.25\n"
        "    &INTERACTION_POTENTIAL\n"
        "      POTENTIAL_TYPE SHORTRANGE\n"
        "      OMEGA 0.11\n"
        "    &END INTERACTION_POTENTIAL\n"
        "    &SCREENING\n"
        "      EPS_SCHWARZ 1.0E-6\n"
        "      SCREEN_ON_INITIAL_P TRUE\n"
        "    &END SCREENING\n"
        "  &END HF\n"
        "&END XC"
    )


def build_admm_block(functional: str) -> str:
    if functional in ("PBE0", "HSE06"):
        return (
            "&AUXILIARY_DENSITY_MATRIX_METHOD\n"
            "  METHOD BASIS_PROJECTION\n"
            "  ADMM_PURIFICATION_METHOD MO_DIAG\n"
            "&END AUXILIARY_DENSITY_MATRIX_METHOD\n"
        )
    return ""


def is_gamma_only(kpoints) -> bool:
    kx, ky, kz = kpoints
    return kx == 1 and ky == 1 and kz == 1


def build_kpoints_block(kpoints) -> str:
    """&KPOINTS block, or nothing for a 1x1x1 mesh.

    Writing &KPOINTS at all (even MONKHORST-PACK 1 1 1) switches CP2K to the
    k-point code path, where &OT aborts ("OT not possible with kpoint
    calculations") -- so a Gamma-only run must omit the section entirely.
    """
    if is_gamma_only(kpoints):
        return ""
    kx, ky, kz = kpoints
    return (
        "    &KPOINTS\n"
        f"      SCHEME MONKHORST-PACK {kx} {ky} {kz}\n"
        "    &END KPOINTS\n"
    )


def build_scf_block(kpoints, disable_wfn_restart: bool,
                     eps_scf: float = 1.0E-6, max_scf: int = 300,
                     outer_max_scf: int = 50, added_mos: int = 0,
                     wfn_restart_file: str = None, mixing_alpha: float = 0.4) -> str:
    """Build the &SCF block, choosing the minimizer CP2K actually supports.

    CP2K's OT (Orbital Transformation) minimizer only supports Gamma-point
    (MO derivatives are unavailable for k-point calculations -- see
    cp2k.org/faq:kpoints); any non-Gamma Monkhorst-Pack mesh MUST use
    traditional DIAGONALIZATION + density MIXING instead, or CP2K aborts
    with "OT not possible with kpoint calculations".

    eps_scf/max_scf/outer_max_scf default to values mirroring a VASP
    production INCAR (EDIFF=1.0E-6, NELM=1500) -- NELM itself is NOT copied
    1:1 since CP2K's OT/DIIS converges far faster than VASP's blocked
    Davidson scheme; max_scf=300 is a generous safety margin, not a forced
    match to 1500.

    added_mos: extra empty bands included in the main SCF diagonalization
    -- needed so the DOS/PDOS actually shows conduction-band states, not
    just occupied valence states. 0 (default) matches normal relaxation
    runs where only occupied states matter.

    wfn_restart_file: if given, seeds SCF_GUESS RESTART from a previously
    converged .wfn (e.g. from the Etot run), instead of SCF_GUESS ATOMIC.
    """
    restart_block = (
        "\n      &PRINT\n        &RESTART OFF\n        &END RESTART\n      &END PRINT"
        if disable_wfn_restart else
        "\n      &PRINT\n        &RESTART ON\n"
        "          BACKUP_COPIES 1\n"
        "        &END RESTART\n      &END PRINT"
    )
    added_mos_line = f"      ADDED_MOS {added_mos}\n" if added_mos else ""
    if wfn_restart_file:
        scf_guess = "RESTART"
        wfn_line = f"      WFN_RESTART_FILE_NAME {wfn_restart_file}\n"
    else:
        scf_guess = "ATOMIC"
        wfn_line = ""

    if is_gamma_only(kpoints):
        return (
            "&SCF\n"
            f"      EPS_SCF {eps_scf}\n"
            f"      MAX_SCF {max_scf}\n"
            f"      SCF_GUESS {scf_guess}\n"
            f"{wfn_line}"
            f"{added_mos_line}"
            "      &OT\n"
            "        MINIMIZER DIIS\n"
            "        PRECONDITIONER FULL_SINGLE_INVERSE\n"
            "      &END OT\n"
            "      &OUTER_SCF\n"
            f"        MAX_SCF {outer_max_scf}\n"
            f"        EPS_SCF {eps_scf}\n"
            f"      &END OUTER_SCF{restart_block}\n"
            "    &END SCF"
        )

    # Non-Gamma k-mesh: OT is unavailable, use DIAGONALIZATION + MIXING.
    # ALPHA/BETA/NBROYDEN below are reasonable starting points, not
    # universally converged -- tune per system if SCF struggles.
    # ALPHA 0.4 mirrors AMIX from the VASP INCAR (same role: density mixing
    # fraction); BETA/NBROYDEN don't have a direct VASP equivalent since the
    # mixing algorithms differ (Broyden here vs VASP's Kerker/RMM-DIIS).
    return (
        "&SCF\n"
        f"      EPS_SCF {eps_scf}\n"
        f"      MAX_SCF {max_scf}\n"
        f"      SCF_GUESS {scf_guess}\n"
        f"{wfn_line}"
        f"{added_mos_line}"
        "      &DIAGONALIZATION\n"
        "        ALGORITHM STANDARD\n"
        "      &END DIAGONALIZATION\n"
        "      &MIXING\n"
        "        METHOD BROYDEN_MIXING\n"
        f"        ALPHA {mixing_alpha}\n"
        "        BETA 1.5\n"
        "        NBROYDEN 8\n"
        "      &END MIXING\n"
        "      &OUTER_SCF\n"
        f"        MAX_SCF {outer_max_scf}\n"
        f"        EPS_SCF {eps_scf}\n"
        f"      &END OUTER_SCF{restart_block}\n"
        "    &END SCF"
    )


def total_valence_electrons(structure, overrides: dict) -> int:
    """Sum GTH valence electrons over all atoms -- the CP2K analog of VASP's
    NELECT, used to auto-size ADDED_MOS the same way NBANDS=NELECT does."""
    total = 0
    for site in structure.sites:
        element = site.specie.symbol
        override = overrides.get(element, {})
        if "valence_electrons" in override:
            total += int(override["valence_electrons"])
            continue
        valence = DEFAULT_GTH_VALENCE.get(element)
        if valence is None:
            sys.exit(
                f"ERROR: no default GTH valence known for element '{element}'. "
                f"Add it to DEFAULT_GTH_VALENCE or supply --kind-overrides with "
                f"a 'valence_electrons' entry for this element."
            )
        total += valence
    return total


def render_input(structure: Structure, args) -> str:
    use_admm = args.functional in ("PBE0", "HSE06")
    overrides = load_kind_overrides(args.kind_overrides)

    # getattr fallbacks: scripts that build a lightweight Namespace directly
    # (e.g. kmesh_stress_convergence.py) may not set these explicitly.
    eps_scf = getattr(args, "eps_scf", 1.0E-6)
    max_scf = getattr(args, "max_scf", 300)
    outer_max_scf = getattr(args, "outer_max_scf", 50)
    max_iter = getattr(args, "max_iter", 500)
    max_force = getattr(args, "max_force", 0.010)
    mixing_alpha = getattr(args, "mixing_alpha", 0.4)
    surface_dipole = getattr(args, "surface_dipole", False)

    cell_block = build_cell_block(structure)
    coord_block = build_coord_block(structure)
    kind_blocks = build_kind_blocks(structure, args.basis_family, args.functional, overrides, use_admm)
    xc_block = build_xc_block(args.functional)
    admm_block = build_admm_block(args.functional)
    motion_block = build_motion_block(args.run_type, max_iter, max_force)

    is_cell_opt = args.run_type == "CELL_OPT"
    # REL_CUTOFF safety margin for the stress tensor (CELL_OPT only) -- see
    # --stress-rel-cutoff-factor help text for the rationale.
    effective_rel_cutoff = (
        round(args.rel_cutoff * args.stress_rel_cutoff_factor) if is_cell_opt else args.rel_cutoff
    )
    stress_tensor_line = "    STRESS_TENSOR ANALYTICAL\n" if is_cell_opt else ""

    run_type_global = "GEO_OPT" if args.run_type in ("GEO_OPT",) else \
                       "CELL_OPT" if args.run_type == "CELL_OPT" else "ENERGY_FORCE"

    is_one_shot = args.run_type in ("ENERGY", "ENERGY_FORCE")
    # One-shot scan points (cutoff/k-mesh single-point convergence tests) never
    # get restarted, so the automatic MO/wavefunction restart write is pure
    # overhead -- and a real liability when many array tasks write it
    # concurrently to a shared filesystem (observed: SIGABRT inside
    # write_mo_set_to_restart under I/O contention). Disable it for those;
    # keep it ON for GEO_OPT/CELL_OPT where it speeds up SCF on resume.
    scf_block = build_scf_block(args.kpoints, disable_wfn_restart=is_one_shot,
                                 eps_scf=eps_scf, max_scf=max_scf, outer_max_scf=outer_max_scf,
                                 mixing_alpha=mixing_alpha)
    kpoints_block = build_kpoints_block(args.kpoints)
    dipole_lines = (
        "    SURFACE_DIPOLE_CORRECTION TRUE\n"
        "    SURF_DIP_DIR Z\n"
    ) if surface_dipole else ""

    template = f"""&GLOBAL
  PROJECT  {args.project_name}
  RUN_TYPE {run_type_global}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
{stress_tensor_line}  &DFT
    BASIS_SET_FILE_NAME  {args.basis_set_file}
    POTENTIAL_FILE_NAME  {args.potential_file}
{dipole_lines}    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
    &END QS
    &MGRID
      CUTOFF {args.cutoff}
      REL_CUTOFF {effective_rel_cutoff}
      NGRIDS 5
    &END MGRID
    {xc_block}
    {scf_block}
{kpoints_block}    {admm_block}&END DFT
  &SUBSYS
    {cell_block}
    {coord_block}
    {kind_blocks}
  &END SUBSYS
&END FORCE_EVAL

{motion_block}"""

    return template


def main():
    args = parse_args()
    structure = get_structure(mp_id=args.mp_id, cif_path=args.cif, api_key=args.api_key)
    input_text = render_input(structure, args)

    output_path = Path(args.output) if args.output else Path.cwd() / f"{args.project_name}.inp"
    output_path.write_text(input_text, encoding="utf-8")

    print(f"Structure: {structure.composition.reduced_formula} "
          f"({len(structure)} sites, functional={args.functional}, run_type={args.run_type})")
    print(f"CP2K input written to: {output_path}")


if __name__ == "__main__":
    main()