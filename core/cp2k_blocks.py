"""
cp2k_blocks.py -- shared CP2K input-block library for the pipeline.

Single source of truth for every block the step scripts (01_..07_) write.
All tunable physics lives in DEFAULTS; step scripts expose it through
add_common_args() and never hardcode values. Only jobs/*.sh carry
cluster-specific values.

Design rules (see CP2K_PIPELINE.md, PART 2):
  * Gamma-only  -> &OT (+ &OUTER_SCF), no &KPOINTS section at all
                   (writing &KPOINTS, even 1 1 1, switches CP2K to the
                   k-point code path, where OT aborts).
  * k != Gamma  -> &DIAGONALIZATION (EPS_ADAPT) + &MIXING + &SMEAR
                   (FERMI_DIRAC, 300 K) + ADDED_MOS.
  * ADDED_MOS is a fixed default (10), NEVER derived from the electron count.
  * &QS EXTRAPOLATION USE_GUESS (robust restarts along GEO/CELL_OPT).
  * &KPOINTS PARALLEL_GROUP_SIZE -1 (one k-point group per MPI rank block).
  * CELL_OPT gets STRESS_TENSOR ANALYTICAL and a &CELL_REF (grids built on a
    larger reference cell, so the effective cutoff never drops as the cell
    breathes).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

try:
    from pymatgen.core import Lattice, Structure
except ImportError:  # pragma: no cover
    sys.exit("ERROR: pymatgen is required (conda activate cp2k_env).")

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903

DEFAULTS = {
    # grids
    "cutoff": 600,            # Ry, MGRID CUTOFF (placeholder until 01 converges it)
    "rel_cutoff": 60,         # Ry, MGRID REL_CUTOFF
    "ngrids": 5,
    "eps_default": 1.0e-12,
    "extrapolation": "USE_GUESS",
    # SCF
    "eps_scf": 1.0e-6,
    "max_scf": 300,           # inner loop (OT or diagonalization)
    "outer_max_scf": 20,      # OT only
    "added_mos": 10,
    "smear_method": "FERMI_DIRAC",
    "smear_temperature": 300.0,  # K
    "eps_adapt": 0.01,
    "mixing_method": "BROYDEN_MIXING",
    "mixing_alpha": 0.4,
    "nbroyden": 8,
    # optimizers
    "optimizer": "BFGS",
    "max_iter": 500,
    "max_force": 0.010,       # eV/angstrom (VASP EDIFFG = -0.01)
    "pressure_tolerance": 100.0,  # bar
    "cell_ref_factor": 1.15,  # linear factor on the lattice vectors for &CELL_REF
    # k-points
    "parallel_group_size": -1,
    # files
    "basis_files": ["BASIS_MOLOPT", "BASIS_MOLOPT_UCL"],
    "potential_file": "GTH_POTENTIALS",
    "potential_family": "GTH-PBE",
}

# GTH valence ("-qN") per element for the MOLOPT/GTH-PBE family shipped with CP2K.
GTH_VALENCE = {
    "H": 1, "Li": 3, "Be": 4, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7,
    "Na": 9, "Mg": 10, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7,
    "K": 9, "Ca": 10, "Sc": 11, "Ti": 12, "V": 13, "Cr": 14, "Mn": 15,
    "Fe": 16, "Co": 17, "Ni": 18, "Cu": 11, "Zn": 12, "Se": 6,
    "Sr": 10, "Y": 11, "Zr": 12, "Nb": 13, "Mo": 14, "Ag": 11,
    "Sn": 4, "Ba": 10, "Hf": 12, "Pb": 4,
}


# ----------------------------------------------------------------------------
# CLI helpers
# ----------------------------------------------------------------------------
def add_common_args(p: argparse.ArgumentParser, *, need_cutoff=True, need_kpoints=True):
    """Arguments every generator shares. --basis is mandatory on purpose:
    the basis family is a physics decision, never a silent default."""
    p.add_argument("--structure", required=True,
                   help="Input structure (CIF, POSCAR, CP2K .restart ...).")
    p.add_argument("--basis", required=True,
                   help="Basis family WITHOUT the -qN suffix, e.g. DZVP-MOLOPT-SR-GTH.")
    p.add_argument("--project", required=True, help="PROJECT name / file stem.")
    p.add_argument("--output-dir", default=".", help="Where to write (default: cwd).")
    p.add_argument("--basis-file", nargs="+", default=DEFAULTS["basis_files"])
    p.add_argument("--potential-file", default=DEFAULTS["potential_file"])
    p.add_argument("--potential-family", default=DEFAULTS["potential_family"])
    if need_cutoff:
        p.add_argument("--cutoff", type=float, default=DEFAULTS["cutoff"])
        p.add_argument("--rel-cutoff", type=float, default=DEFAULTS["rel_cutoff"])
    p.add_argument("--ngrids", type=int, default=DEFAULTS["ngrids"])
    p.add_argument("--eps-default", type=float, default=DEFAULTS["eps_default"])
    if need_kpoints:
        p.add_argument("--kpoints", type=int, nargs=3, default=[1, 1, 1],
                       metavar=("KX", "KY", "KZ"))
    p.add_argument("--eps-scf", type=float, default=DEFAULTS["eps_scf"])
    p.add_argument("--max-scf", type=int, default=DEFAULTS["max_scf"])
    p.add_argument("--outer-max-scf", type=int, default=DEFAULTS["outer_max_scf"])
    p.add_argument("--added-mos", type=int, default=DEFAULTS["added_mos"])
    p.add_argument("--smear-temperature", type=float, default=DEFAULTS["smear_temperature"])
    p.add_argument("--mixing-alpha", type=float, default=DEFAULTS["mixing_alpha"])
    p.add_argument("--eps-adapt", type=float, default=DEFAULTS["eps_adapt"])
    return p


# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------
def read_structure(path: str | Path) -> Structure:
    path = Path(path)
    if not path.exists():
        sys.exit(f"ERROR: structure file not found: {path}")
    if path.suffix == ".restart" or path.name.endswith(".restart") or ".restart.bak" in path.name:
        return read_restart_structure(path)
    return Structure.from_file(str(path))


def read_restart_structure(path: str | Path) -> Structure:
    """Final cell + coordinates from a CP2K .restart (the CONTCAR analog).

    Handles A/B/C or ABC(+ALPHA_BETA_GAMMA) cells, SCALED or Cartesian
    (angstrom) &COORD, and '&KIND Zr' element labels with suffixes (Zr_1).
    """
    text = Path(path).read_text()
    subsys = text[text.upper().find("&SUBSYS"):]
    cell_txt = re.search(r"&CELL\b(.*?)&END CELL", subsys, re.S | re.I).group(1)
    cell_txt = re.sub(r"&CELL_REF.*?&END CELL_REF", "", cell_txt, flags=re.S | re.I)
    vecs = {}
    for key in ("A", "B", "C"):
        m = re.search(rf"^\s*{key}\s+(?:\[[^\]]*\]\s+)?([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)",
                      cell_txt, re.M)
        if m:
            vecs[key] = [float(x) for x in m.groups()]
    if len(vecs) == 3:
        lattice = Lattice([vecs["A"], vecs["B"], vecs["C"]])
    else:
        abc = re.search(r"^\s*ABC\s+(?:\[[^\]]*\]\s+)?([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)",
                        cell_txt, re.M)
        ang = re.search(r"^\s*ALPHA_BETA_GAMMA\s+(?:\[[^\]]*\]\s+)?([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)",
                        cell_txt, re.M)
        a, b, c = (float(x) for x in abc.groups())
        al, be, ga = (float(x) for x in ang.groups()) if ang else (90.0, 90.0, 90.0)
        lattice = Lattice.from_parameters(a, b, c, al, be, ga)

    coord_txt = re.search(r"&COORD\b(.*?)&END COORD", subsys, re.S | re.I).group(1)
    scaled = bool(re.search(r"^\s*SCALED\s*(\.TRUE\.|T|TRUE)?\s*$", coord_txt, re.M | re.I))
    species, coords = [], []
    for line in coord_txt.splitlines():
        parts = line.split()
        if len(parts) < 4 or parts[0].upper() in ("SCALED", "UNIT"):
            continue
        try:
            xyz = [float(x) for x in parts[1:4]]
        except ValueError:
            continue
        species.append(re.match(r"[A-Z][a-z]?", parts[0]).group(0))
        coords.append(xyz)
    return Structure(lattice, species, coords, coords_are_cartesian=not scaled)


def n_formula_units(structure: Structure) -> int:
    return int(round(structure.composition.get_reduced_composition_and_factor()[1]))


def kmesh_from_density(lattice: Lattice, density: float, periodic=(True, True, True)):
    """k_i = ceil(density / |a_i|) -- density is n_i*|a_i| in angstrom."""
    return [max(1, math.ceil(density / L)) if per else 1
            for L, per in zip(lattice.abc, periodic)]


def kmesh_density(lattice: Lattice, kpoints) -> float:
    """Largest n_i*|a_i| of a mesh (the value to reuse on other cells)."""
    return max(k * L for k, L in zip(kpoints, lattice.abc))


# ----------------------------------------------------------------------------
# Blocks
# ----------------------------------------------------------------------------
def is_gamma_only(kpoints) -> bool:
    return list(kpoints) == [1, 1, 1]


def build_global_block(project: str, run_type: str, print_level: str = "MEDIUM") -> str:
    return (f"&GLOBAL\n  PROJECT {project}\n  RUN_TYPE {run_type}\n"
            f"  PRINT_LEVEL {print_level}\n&END GLOBAL\n")


def build_cell_block(structure: Structure, cell_ref_factor: float | None = None,
                     periodic: str = "XYZ") -> str:
    """&CELL with explicit A/B/C vectors. cell_ref_factor (CELL_OPT only) adds
    a &CELL_REF with every lattice vector scaled by that linear factor."""
    m = structure.lattice.matrix
    lines = ["    &CELL"]
    for lab, v in zip("ABC", m):
        lines.append(f"      {lab} [angstrom] {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
    lines.append(f"      PERIODIC {periodic}")
    if cell_ref_factor:
        lines.append("      &CELL_REF")
        for lab, v in zip("ABC", m * cell_ref_factor):
            lines.append(f"        {lab} [angstrom] {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
        lines.append("      &END CELL_REF")
    lines.append("    &END CELL")
    return "\n".join(lines) + "\n"


def build_coord_block(structure: Structure) -> str:
    lines = ["    &COORD", "      SCALED"]
    for site in structure:
        fx, fy, fz = site.frac_coords
        lines.append(f"      {site.specie.symbol:<3} {fx: .10f} {fy: .10f} {fz: .10f}")
    lines.append("    &END COORD")
    return "\n".join(lines) + "\n"


def basis_name(element: str, basis: str) -> str:
    q = GTH_VALENCE.get(element)
    if q is None:
        sys.exit(f"ERROR: no GTH valence known for {element}; add it to GTH_VALENCE.")
    return f"{basis}-q{q}"


def build_kind_blocks(structure: Structure, basis: str, potential_family: str = "GTH-PBE") -> str:
    out = []
    for el in sorted({s.specie.symbol for s in structure}):
        bname = basis_name(el, basis)  # exits if the element has no GTH valence
        out.append(f"    &KIND {el}\n      BASIS_SET {bname}\n"
                   f"      POTENTIAL {potential_family}-q{GTH_VALENCE[el]}\n    &END KIND")
    return "\n".join(out) + "\n"


def valence_electrons(structure: Structure) -> int:
    return sum(GTH_VALENCE[s.specie.symbol] for s in structure)


def build_qs_block(eps_default: float = DEFAULTS["eps_default"],
                   extrapolation: str = DEFAULTS["extrapolation"]) -> str:
    return (f"    &QS\n      METHOD GPW\n      EPS_DEFAULT {eps_default:.1E}\n"
            f"      EXTRAPOLATION {extrapolation}\n    &END QS\n")


def build_mgrid_block(cutoff: float, rel_cutoff: float, ngrids: int = DEFAULTS["ngrids"]) -> str:
    return (f"    &MGRID\n      CUTOFF {cutoff:g}\n      REL_CUTOFF {rel_cutoff:g}\n"
            f"      NGRIDS {ngrids}\n    &END MGRID\n")


def build_xc_block(functional: str = "PBE") -> str:
    if functional.upper() != "PBE":
        sys.exit("ERROR: only PBE is implemented in this pipeline copy "
                 "(hybrid/RI-HFXk work lives on coaraci, see SECOND_BRAIN.md).")
    return "    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n"


def build_scf_block(kpoints, *, eps_scf=DEFAULTS["eps_scf"], max_scf=DEFAULTS["max_scf"],
                    outer_max_scf=DEFAULTS["outer_max_scf"], added_mos=DEFAULTS["added_mos"],
                    smear_temperature=DEFAULTS["smear_temperature"],
                    mixing_alpha=DEFAULTS["mixing_alpha"], eps_adapt=DEFAULTS["eps_adapt"],
                    scf_guess="ATOMIC", wfn_restart=None, restart_print=True,
                    outer_scf=True, ignore_convergence_failure=False) -> str:
    """OT for Gamma-only, DIAGONALIZATION+MIXING+SMEAR for k-points.

    restart_print: &SCF%PRINT%RESTART ON/OFF, always written explicitly
    (OFF for concurrent array scans: the binary k-point restart writer
    crashed under shared-FS I/O contention in the legacy pipeline).
    outer_scf: False for the grid scan (MAX_SCF 1, no outer loop).
    ignore_convergence_failure: CP2K 2026.1 ABORTS on an unconverged SCF by
    default; only the grid scan (deliberately MAX_SCF 1) sets this.
    """
    L = ["    &SCF", f"      EPS_SCF {eps_scf:.1E}", f"      MAX_SCF {max_scf}",
         f"      SCF_GUESS {scf_guess}"]
    if ignore_convergence_failure:
        L.append("      IGNORE_CONVERGENCE_FAILURE TRUE")
    if wfn_restart:
        L.append(f"      WFN_RESTART_FILE_NAME {wfn_restart}")
    if is_gamma_only(kpoints):
        L += ["      &OT", "        MINIMIZER DIIS", "        PRECONDITIONER FULL_SINGLE_INVERSE",
              "      &END OT"]
        if outer_scf:
            L += ["      &OUTER_SCF", f"        MAX_SCF {outer_max_scf}",
                  f"        EPS_SCF {eps_scf:.1E}", "      &END OUTER_SCF"]
    else:
        L += [f"      ADDED_MOS {added_mos}",
              "      &DIAGONALIZATION", "        ALGORITHM STANDARD",
              f"        EPS_ADAPT {eps_adapt}", "      &END DIAGONALIZATION",
              "      &MIXING", f"        METHOD {DEFAULTS['mixing_method']}",
              f"        ALPHA {mixing_alpha}", f"        NBROYDEN {DEFAULTS['nbroyden']}",
              "      &END MIXING",
              "      &SMEAR ON", f"        METHOD {DEFAULTS['smear_method']}",
              f"        ELECTRONIC_TEMPERATURE [K] {smear_temperature:g}", "      &END SMEAR"]
    L += ["      &PRINT", f"        &RESTART {'ON' if restart_print else 'OFF'}",
          "          BACKUP_COPIES 1" if restart_print else None, "        &END RESTART",
          "      &END PRINT", "    &END SCF"]
    return "\n".join(x for x in L if x is not None) + "\n"


def build_kpoints_block(kpoints, parallel_group_size=DEFAULTS["parallel_group_size"]) -> str:
    if is_gamma_only(kpoints):
        return ""
    kx, ky, kz = kpoints
    return (f"    &KPOINTS\n      SCHEME MONKHORST-PACK {kx} {ky} {kz}\n"
            f"      PARALLEL_GROUP_SIZE {parallel_group_size}\n    &END KPOINTS\n")


def build_motion_block(run_type: str, *, max_iter=DEFAULTS["max_iter"],
                       max_force=DEFAULTS["max_force"], optimizer=DEFAULTS["optimizer"],
                       pressure_tolerance=DEFAULTS["pressure_tolerance"]) -> str:
    if run_type not in ("GEO_OPT", "CELL_OPT"):
        return ""
    sec = run_type
    body = [f"  &{sec}"]
    if run_type == "CELL_OPT":
        body += ["    TYPE DIRECT_CELL_OPT", f"    PRESSURE_TOLERANCE [bar] {pressure_tolerance:g}"]
    body += [f"    OPTIMIZER {optimizer}", f"    MAX_ITER {max_iter}",
             f"    MAX_FORCE [eV*angstrom^-1] {max_force}", f"  &END {sec}"]
    return ("&MOTION\n" + "\n".join(body) + "\n"
            "  &PRINT\n    &RESTART\n      BACKUP_COPIES 3\n      &EACH\n"
            f"        {sec} 1\n      &END EACH\n    &END RESTART\n  &END PRINT\n&END MOTION\n")


def build_input(structure: Structure, args, *, run_type: str, kpoints,
                cell_ref_factor=None, stress=False, surface_dipole=False,
                dft_print="", print_level="MEDIUM", scf_kwargs=None,
                motion_kwargs=None) -> str:
    """Assemble the full input. dft_print: text of an &PRINT section for &DFT."""
    scf_kwargs = dict(scf_kwargs or {})
    scf = build_scf_block(
        kpoints, eps_scf=args.eps_scf, max_scf=scf_kwargs.pop("max_scf", args.max_scf),
        outer_max_scf=args.outer_max_scf, added_mos=args.added_mos,
        smear_temperature=args.smear_temperature, mixing_alpha=args.mixing_alpha,
        eps_adapt=args.eps_adapt, **scf_kwargs)
    basis_lines = "".join(f"    BASIS_SET_FILE_NAME {b}\n" for b in args.basis_file)
    dip = "    SURFACE_DIPOLE_CORRECTION TRUE\n    SURF_DIP_DIR Z\n" if surface_dipole else ""
    return (
        build_global_block(args.project, run_type, print_level) + "\n"
        + "&FORCE_EVAL\n  METHOD Quickstep\n"
        + ("  STRESS_TENSOR ANALYTICAL\n" if stress else "")
        + "  &DFT\n" + basis_lines
        + f"    POTENTIAL_FILE_NAME {args.potential_file}\n" + dip
        + build_qs_block(args.eps_default)
        + build_mgrid_block(args.cutoff, args.rel_cutoff, args.ngrids)
        + build_xc_block("PBE") + scf + build_kpoints_block(kpoints) + dft_print
        + "  &END DFT\n  &SUBSYS\n"
        + build_cell_block(structure, cell_ref_factor)
        + build_coord_block(structure)
        + build_kind_blocks(structure, args.basis, args.potential_family)
        + "  &END SUBSYS\n&END FORCE_EVAL\n\n"
        + build_motion_block(run_type, **(motion_kwargs or {}))
    )


def write_input(text: str, outdir: str | Path, project: str) -> Path:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{project}.inp"
    path.write_text(text)
    return path
