#!/usr/bin/env python3
"""
make_benchmark.py -- CP2K functional benchmark for Ag2SeO3 (two-temperature model, VASP -> CP2K).

Reproduces the student's VASP setup (INCAR_DOS: ISMEAR=0 Gaussian smearing, SIGMA scan,
NBANDS=300, EDIFF=1E-6, EDIFFG=-0.01) with four GGA functionals: PBE (reference),
PBEsol, revPBE, BLYP.

Stages (same logic as the cp2k-pipeline 2D method; smearing fixed at SIGMA=0.01 eV there):
  kmesh           k-mesh convergence, CELL_OPT per mesh (PBE, sigma 0.01), after grid
  grid-cutoff     CUTOFF scan at fixed REL_CUTOFF         (ENERGY, Gamma, MAX_SCF 1)
  grid-relcutoff  REL_CUTOFF scan at the chosen CUTOFF    (ENERGY, Gamma, MAX_SCF 1)
  etot            "Etot" = GEO_OPT with FIXED cell (VASP ISIF=2) for every SIGMA

Only the Etot stage varies SIGMA. Smearing is GAUSSIAN (ISMEAR=0) in every stage.

Needs only numpy-free pure Python (no pymatgen on harpia).

Usage:
  python make_benchmark.py grid-cutoff    --poscar vasp_inputs/POSCAR_Convencional --out-dir .
  python make_benchmark.py grid-relcutoff --poscar ... --cutoff PBE=600 PBEsol=600 revPBE=600 BLYP=600
  python make_benchmark.py etot           --poscar ... --cutoff PBE=600 ... --rel-cutoff PBE=50 ...
"""
import argparse
import sys
from pathlib import Path

# functional -> (GTH potential family, XC_FUNCTIONAL body)
FUNCTIONALS = {
    "PBE": ("GTH-PBE", "      &PBE\n        PARAMETRIZATION ORIG\n      &END PBE\n"),
    "PBEsol": ("GTH-PBE", "      &PBE\n        PARAMETRIZATION PBESOL\n      &END PBE\n"),
    "revPBE": ("GTH-PBE", "      &PBE\n        PARAMETRIZATION REVPBE\n      &END PBE\n"),
    "BLYP": ("GTH-BLYP", "      &BECKE88\n      &END BECKE88\n      &LYP\n      &END LYP\n"),
}
GTH_Q = {"Ag": 11, "Se": 6, "O": 6}          # valence electrons of the GTH potentials
BASIS = "DZVP-MOLOPT-SR-GTH"                  # exists for Ag, Se and O (BASIS_MOLOPT)
SIGMAS = [0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 6]   # eV, student's scan
NBANDS_VASP = 300                             # INCAR NBANDS


def read_poscar(path):
    L = [l.rstrip("\n") for l in Path(path).read_text().splitlines()]
    scale = float(L[1].split()[0])
    lat = [[float(x) * scale for x in L[i].split()[:3]] for i in (2, 3, 4)]
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    assert L[7].strip().lower().startswith("d"), "expected Direct coordinates"
    n = sum(cnt)
    frac = [[float(x) for x in L[8 + i].split()[:3]] for i in range(n)]
    species = [e for e, c in zip(els, cnt) for _ in range(c)]
    return lat, species, frac


def added_mos(sigma, nelec, nbands=NBANDS_VASP):
    """ADDED_MOS so that the band count tracks VASP NBANDS=300 where the smearing needs it.
    nocc = nelec/2. sigma<=0.1 eV: few empty bands; sigma>=1 eV: full VASP NBANDS."""
    nocc = nelec // 2
    full = nbands - nocc
    f = max(1, round(nelec / 184))          # 184 e- = 24-atom conventional cell
    if sigma <= 0.1:
        return 20 * f
    if sigma <= 0.5:
        return 60 * f
    if sigma <= 1:
        return 120 * f
    return full


def block_input(project, run_type, func, lat, species, frac, cutoff, rel_cutoff, sigma,
                kpts, grid_scan=False, eps_scf=1.0e-7, nbands=NBANDS_VASP, md=None):
    cell_opt = run_type == "CELL_OPT"
    pot, xc = FUNCTIONALS[func]
    nelec = sum(GTH_Q[s] for s in species)
    mos = 20 if grid_scan else added_mos(sigma, nelec, nbands)
    gamma = tuple(kpts) == (1, 1, 1) or grid_scan
    scf = ["    &SCF", f"      EPS_SCF {eps_scf:.1E}",
           "      MAX_SCF 1" if grid_scan else "      MAX_SCF 1500", "      SCF_GUESS ATOMIC"]
    if grid_scan:
        scf.append("      IGNORE_CONVERGENCE_FAILURE TRUE")
    scf += [f"      ADDED_MOS {mos}",
            "      &DIAGONALIZATION", "        ALGORITHM STANDARD", "        EPS_ADAPT 0.01",
            "      &END DIAGONALIZATION",
            "      &MIXING", "        METHOD BROYDEN_MIXING", "        ALPHA 0.4", "        NBROYDEN 8",
            "      &END MIXING",
            "      &SMEAR ON", "        METHOD GAUSSIAN", f"        SIGMA [eV] {sigma:g}", "      &END SMEAR",
            "      &PRINT", "        &RESTART OFF", "        &END RESTART", "      &END PRINT",
            "    &END SCF"]
    kp = "" if gamma else (f"    &KPOINTS\n      SCHEME MONKHORST-PACK {kpts[0]} {kpts[1]} {kpts[2]}\n"
                           "      PARALLEL_GROUP_SIZE -1\n    &END KPOINTS\n")
    cell = "    &CELL\n" + "".join(
        f"      {a} [angstrom] {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n" for a, v in zip("ABC", lat)
    ) + "      PERIODIC XYZ\n" + (
        "    &CELL_REF\n" + "".join(
            f"      {a} [angstrom] {1.15 * v[0]:.10f} {1.15 * v[1]:.10f} {1.15 * v[2]:.10f}\n"
            for a, v in zip("ABC", lat)) + "    &END CELL_REF\n" if cell_opt else "") + "    &END CELL\n"
    coord = "    &COORD\n      SCALED\n" + "".join(
        f"      {s:<3} {f[0]: .10f} {f[1]: .10f} {f[2]: .10f}\n" for s, f in zip(species, frac)
    ) + "    &END COORD\n"
    kinds = "".join(
        f"    &KIND {e}\n      BASIS_SET {BASIS}-q{GTH_Q[e]}\n      POTENTIAL {pot}-q{GTH_Q[e]}\n    &END KIND\n"
        for e in sorted(set(species)))
    out = (f"&GLOBAL\n  PROJECT {project}\n  RUN_TYPE {run_type}\n"
           f"  PRINT_LEVEL {'MEDIUM' if grid_scan else 'LOW'}\n" + (f"  SEED {md['seed']}\n" if md else "") + "&END GLOBAL\n\n"
           "&FORCE_EVAL\n  METHOD Quickstep\n" + ("  STRESS_TENSOR ANALYTICAL\n" if cell_opt else "") + "  &DFT\n"
           "    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n"
           "    &QS\n      METHOD GPW\n      EPS_DEFAULT 1.0E-12\n"
           + (f"      EXTRAPOLATION {'ASPC' if md else 'USE_GUESS'}\n"
              + ("      EXTRAPOLATION_ORDER 3\n" if md else "")) + "    &END QS\n"
           f"    &MGRID\n      CUTOFF {cutoff:g}\n      REL_CUTOFF {rel_cutoff:g}\n      NGRIDS 5\n    &END MGRID\n"
           f"    &XC\n      &XC_FUNCTIONAL\n{xc}      &END XC_FUNCTIONAL\n    &END XC\n"
           + "\n".join(scf) + "\n" + kp + "  &END DFT\n  &SUBSYS\n" + cell + coord + kinds
           + "  &END SUBSYS\n&END FORCE_EVAL\n")
    if run_type == "MD":
        out += (f"\n&MOTION\n  &MD\n    ENSEMBLE NVE\n    STEPS {md['steps']}\n    TIMESTEP [fs] {md['dt']:g}\n"
                f"    TEMPERATURE [K] {md['temp']:g}\n  &END MD\n"
                "  &PRINT\n    &TRAJECTORY\n      &EACH\n        MD 5\n      &END EACH\n    &END TRAJECTORY\n"
                "    &RESTART\n      BACKUP_COPIES 1\n      &EACH\n        MD 50\n      &END EACH\n    &END RESTART\n"
                "  &END PRINT\n&END MOTION\n")
    if run_type == "CELL_OPT":
        out += ("\n&MOTION\n  &CELL_OPT\n    EXTERNAL_PRESSURE [bar] 0\n"
                "    PRESSURE_TOLERANCE [bar] 100\n    OPTIMIZER LBFGS\n    MAX_ITER 800\n"
                "    MAX_FORCE [eV*angstrom^-1] 0.01\n  &END CELL_OPT\n"
                "  &PRINT\n    &RESTART\n      BACKUP_COPIES 1\n      &EACH\n        CELL_OPT 1\n      &END EACH\n"
                "    &END RESTART\n  &END PRINT\n&END MOTION\n")
    if run_type == "GEO_OPT":
        out += ("\n&MOTION\n  &GEO_OPT\n    OPTIMIZER LBFGS\n    MAX_ITER 800\n"
                "    MAX_FORCE [eV*angstrom^-1] 0.01\n  &END GEO_OPT\n"
                "  &PRINT\n    &TRAJECTORY\n      &EACH\n        GEO_OPT 1\n      &END EACH\n    &END TRAJECTORY\n"
                "    &RESTART\n      BACKUP_COPIES 1\n      &EACH\n        GEO_OPT 1\n      &END EACH\n"
                "    &END RESTART\n  &END PRINT\n&END MOTION\n")
    return out


def norm(v):
    return sum(x * x for x in v) ** 0.5


def parse_kv(items, label):
    d = {}
    for it in items or []:
        k, v = it.split("=")
        if k not in FUNCTIONALS:
            sys.exit(f"ERROR: unknown functional {k} in --{label}")
        d[k] = float(v)
    return d


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stage", choices=["grid-cutoff", "grid-relcutoff", "kmesh", "etot", "md"])
    p.add_argument("--poscar", required=True)
    p.add_argument("--out-dir", default=".")
    p.add_argument("--functionals", nargs="+", default=list(FUNCTIONALS), choices=list(FUNCTIONALS))
    p.add_argument("--cutoff", nargs="+", metavar="FUNC=Ry",
                   help="per-functional CUTOFF (grid-relcutoff, etot)")
    p.add_argument("--rel-cutoff", nargs="+", metavar="FUNC=Ry", help="per-functional REL_CUTOFF (etot)")
    p.add_argument("--cutoff-values", type=float, nargs="+",
                   default=[300, 400, 500, 600, 700, 800, 900, 1000])
    p.add_argument("--relcutoff-values", type=float, nargs="+", default=[20, 30, 40, 50, 60, 80, 100])
    p.add_argument("--fixed-rel-cutoff", type=float, default=60.0, help="REL_CUTOFF during CUTOFF scan")
    p.add_argument("--densities", type=float, nargs="+", default=[15, 20, 25, 30, 40],
                   help="kmesh stage: k-density n_i*|a_i| in Angstrom, mesh n_i=ceil(d/|a_i|)")
    p.add_argument("--kmesh-func", default="PBE", choices=list(FUNCTIONALS))
    p.add_argument("--md-steps", type=int, default=10000, help="INCAR NSW=10000")
    p.add_argument("--md-dt", type=float, default=0.5, help="fs, INCAR POTIM=0.5")
    p.add_argument("--md-temp", type=float, default=100.0, help="K, INCAR TEBEG=100")
    p.add_argument("--md-seed", type=int, default=12345, help="same seed -> same initial velocities everywhere")
    p.add_argument("--md-nbands", type=int, default=1600, help="INCAR_MD NBANDS=1600")
    p.add_argument("--sigmas", type=float, nargs="+", default=SIGMAS)
    p.add_argument("--kpoints", type=int, nargs=3, default=[4, 2, 3],  # md: pass 1 1 1
                   help="Monkhorst-Pack mesh for Etot (ASSUMED; VASP KPOINTS not provided)")
    a = p.parse_args()

    lat, sp, fr = read_poscar(a.poscar)
    cut, rel = parse_kv(a.cutoff, "cutoff"), parse_kv(a.rel_cutoff, "rel-cutoff")
    root = Path(a.out_dir)
    runs = []
    for f in a.functionals:
        if a.stage == "grid-cutoff":
            for v in a.cutoff_values:
                runs.append((root / f / "grid" / f"cutoff{v:g}", f"Ag2SeO3_{f}_cutoff{v:g}",
                             "ENERGY", f, v, a.fixed_rel_cutoff, 0.01, True))
        elif a.stage == "grid-relcutoff":
            if f not in cut:
                sys.exit(f"ERROR: --cutoff {f}=Ry required")
            for v in a.relcutoff_values:
                runs.append((root / f / "grid" / f"relcutoff{v:g}", f"Ag2SeO3_{f}_relcutoff{v:g}",
                             "ENERGY", f, cut[f], v, 0.01, True))
        elif a.stage == "md":
            if f not in cut or f not in rel:
                sys.exit(f"ERROR: --cutoff {f}=Ry and --rel-cutoff {f}=Ry required for md")
            for s in a.sigmas:
                runs.append((root / f / "md" / f"sigma{s:g}", f"Ag2SeO3_{f}_md_sigma{s:g}",
                             "MD", f, cut[f], rel[f], s, False, tuple(a.kpoints)))
        elif a.stage == "kmesh":
            if f != a.kmesh_func:
                continue
            if f not in cut or f not in rel:
                sys.exit(f"ERROR: --cutoff {f}=Ry and --rel-cutoff {f}=Ry required for kmesh")
            seen = set()
            for d in a.densities:
                m = tuple(max(1, -(-int(d * 1000) // int(norm(v) * 1000))) for v in lat)
                if m in seen:
                    continue
                seen.add(m)
                runs.append((root / f / "kmesh" / ("k%d_%d_%d" % m), "Ag2SeO3_%s_k%d%d%d" % ((f,) + m),
                             "CELL_OPT", f, cut[f], rel[f], 0.01, False, m))
        else:
            if f not in cut or f not in rel:
                sys.exit(f"ERROR: --cutoff {f}=Ry and --rel-cutoff {f}=Ry required for etot")
            for s in a.sigmas:
                runs.append((root / f / "etot" / f"sigma{s:g}", f"Ag2SeO3_{f}_sigma{s:g}",
                             "GEO_OPT", f, cut[f], rel[f], s, False))
    lst = []
    for run in runs:
        d, proj, rt, f, c, r, s, g = run[:8]
        kp = run[8] if len(run) > 8 else a.kpoints
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{proj}.inp").write_text(block_input(
            proj, rt, f, lat, sp, fr, c, r, s, kp, g,
            eps_scf=4.0e-7 if rt == "MD" else 1.0e-7,
            nbands=a.md_nbands if rt == "MD" else NBANDS_VASP,
            md=dict(steps=a.md_steps, dt=a.md_dt, temp=a.md_temp, seed=a.md_seed) if rt == "MD" else None))
        lst.append(f"{d.resolve()} {proj}.inp")
    listfile = root / f"runs_{a.stage}.txt"
    listfile.write_text("\n".join(lst) + "\n")
    print(f"{len(runs)} inputs; list: {listfile} ({len(runs)} lines)")


if __name__ == "__main__":
    main()
