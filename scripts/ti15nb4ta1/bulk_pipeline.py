#!/usr/bin/env python3
"""
bulk_pipeline.py -- bulk convergence of Nb4Ta1Ti15 (BCC, 20 atoms) in CP2K, one job, all stages.

Pure Python (no numpy/pymatgen: harpia has none). Criteria carried over from the VASP INCAR:
  EDIFF 1E-6 eV -> EPS_SCF 3.67E-8 Ha     NELM 1500 -> MAX_SCF 1500
  EDIFFG -0.01  -> MAX_FORCE 0.01 eV/A    NSW 800   -> MAX_ITER 800
  ISMEAR 0, SIGMA 0.01 -> SMEAR GAUSSIAN SIGMA [eV] 0.01     ISPIN 1 (RKS)    PBE
  KPOINTS 5x4x4 Gamma-centred -> included as a candidate mesh (CP2K MONKHORST-PACK with
  even n is not Gamma-centred; odd n is -- see notes in the report).
Stages (each restart-safe: rerun the same command and finished points are skipped):
  cutoff     ENERGY, Gamma, MAX_SCF 1, CUTOFF scan at fixed REL_CUTOFF   (CP2K 2D method)
  relcutoff  same, REL_CUTOFF scan at the CUTOFF chosen above
  kmesh      CELL_OPT (ions+cell, P=0) per k-mesh -> dE meV/atom, dV %, residual pressure
  cellopt    production CELL_OPT with the chosen CUTOFF/REL_CUTOFF/k-mesh -> bulk.json + POSCAR
Points of a stage run concurrently, each pinned (taskset) to its own slice of the allocated CPUs.

Usage:  python3 bulk_pipeline.py --stage all --ncpus 48
"""
import argparse, glob, json, math, os, re, subprocess, sys, threading, time
from pathlib import Path

HA_EV = 27.211386245988
HERE = Path(__file__).resolve().parent
GTH_Q = {"Ti": 12, "Nb": 13, "Ta": 13}
BASIS = "DZVP-MOLOPT-SR-GTH"
SIGMA_EV = 0.01          # INCAR SIGMA
EPS_SCF = 3.67e-8        # INCAR EDIFF 1e-6 eV in Ha
MAX_SCF = 1500           # INCAR NELM / rule 9
MAX_ITER = 800           # INCAR NSW / rule 9
MAX_FORCE = 0.010        # eV/A, INCAR EDIFFG
ADDED_MOS = 20           # metal, sigma 0.01 eV: a few empty bands (VASP default NBANDS ~ nocc+10)
MESHES = [(2, 2, 2), (3, 3, 3), (4, 3, 3), (5, 4, 4), (6, 5, 5), (7, 6, 6)]  # 5x4x4 = VASP KPOINTS


# ----------------------------------------------------------------- structure (VASP CONTCAR/POSCAR)
def read_poscar(path):
    L = Path(path).read_text().splitlines()
    scale = float(L[1].split()[0])
    lat = [[float(x) * scale for x in L[i].split()[:3]] for i in (2, 3, 4)]
    els, cnt = L[5].split(), [int(x) for x in L[6].split()]
    assert L[7].strip().lower().startswith("d"), "expected Direct coordinates"
    n = sum(cnt)
    frac = [[float(x) for x in L[8 + i].split()[:3]] for i in range(n)]
    sp = [e for e, c in zip(els, cnt) for _ in range(c)]
    return lat, sp, frac


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def write_poscar(path, lat, sp, frac, comment):
    order = []
    for e in sp:
        if e not in order:
            order.append(e)
    with open(path, "w") as f:
        f.write(comment + "\n1.0\n")
        for v in lat:
            f.write("  %.10f %.10f %.10f\n" % tuple(v))
        f.write(" ".join(order) + "\n" + " ".join(str(sp.count(e)) for e in order) + "\nDirect\n")
        for e in order:
            for s, x in zip(sp, frac):
                if s == e:
                    f.write("  %.10f %.10f %.10f\n" % tuple(x))


# ----------------------------------------------------------------- input generation
def make_input(project, run_type, lat, sp, frac, cutoff, rel, kmesh=None, cell_ref=None):
    gamma = kmesh is None
    L = [f"&GLOBAL\n  PROJECT {project}\n  RUN_TYPE {run_type}\n  PRINT_LEVEL MEDIUM\n&END GLOBAL\n",
         "&FORCE_EVAL\n  METHOD Quickstep"]
    if run_type == "CELL_OPT":
        L.append("  STRESS_TENSOR ANALYTICAL")
    scf_extra = ""
    max_scf = MAX_SCF
    if run_type == "ENERGY":      # grid scan: one SCF cycle from the same ATOMIC guess, error cancels
        max_scf = 1
        scf_extra = "      IGNORE_CONVERGENCE_FAILURE TRUE\n"
    L.append("  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n    POTENTIAL_FILE_NAME GTH_POTENTIALS\n"
             "    &QS\n      METHOD GPW\n      EPS_DEFAULT 1.0E-12\n      EXTRAPOLATION USE_GUESS\n    &END QS\n"
             f"    &MGRID\n      CUTOFF {cutoff:g}\n      REL_CUTOFF {rel:g}\n      NGRIDS 5\n    &END MGRID\n"
             "    &XC\n      &XC_FUNCTIONAL PBE\n      &END XC_FUNCTIONAL\n    &END XC\n"
             # 245 valence electrons (odd) + metal + Gaussian smearing => diagonalization, never OT
             f"    &SCF\n      EPS_SCF {EPS_SCF:.2E}\n      MAX_SCF {max_scf}\n      SCF_GUESS ATOMIC\n{scf_extra}"
             f"      ADDED_MOS {ADDED_MOS}\n"
             "      &DIAGONALIZATION\n        ALGORITHM STANDARD\n        EPS_ADAPT 0.01\n      &END DIAGONALIZATION\n"
             "      &MIXING\n        METHOD BROYDEN_MIXING\n        ALPHA 0.4\n        NBROYDEN 8\n      &END MIXING\n"
             f"      &SMEAR ON\n        METHOD GAUSSIAN\n        SIGMA [eV] {SIGMA_EV:g}\n      &END SMEAR\n"
             "      &PRINT\n        &RESTART OFF\n        &END RESTART\n      &END PRINT\n    &END SCF")
    if not gamma:
        L.append("    &KPOINTS\n      SCHEME MONKHORST-PACK %d %d %d\n      PARALLEL_GROUP_SIZE -1\n    &END KPOINTS" % tuple(kmesh))
    L.append("  &END DFT\n  &SUBSYS\n    &CELL")
    for lab, v in zip("ABC", lat):
        L.append("      %s [angstrom] %.10f %.10f %.10f" % ((lab,) + tuple(v)))
    L.append("      PERIODIC XYZ")
    if cell_ref:
        L.append("      &CELL_REF")
        for lab, v in zip("ABC", lat):
            L.append("        %s [angstrom] %.10f %.10f %.10f" % ((lab,) + tuple(x * cell_ref for x in v)))
        L.append("      &END CELL_REF")
    L.append("    &END CELL\n    &COORD\n      SCALED")
    for s, x in zip(sp, frac):
        L.append("      %-3s %.10f %.10f %.10f" % ((s,) + tuple(x)))
    L.append("    &END COORD")
    for el in sorted(set(sp)):
        L.append(f"    &KIND {el}\n      BASIS_SET {BASIS}-q{GTH_Q[el]}\n      POTENTIAL GTH-PBE-q{GTH_Q[el]}\n    &END KIND")
    L.append("  &END SUBSYS\n&END FORCE_EVAL\n")
    if run_type == "CELL_OPT":
        # 2026.2 removed CELL_OPT/TYPE (always DIRECT_CELL_OPT): do not write it.
        # EXTERNAL_PRESSURE explicit: the CP2K default is 100 bar. 20 atoms > 10 -> LBFGS.
        L.append("&MOTION\n  &CELL_OPT\n    EXTERNAL_PRESSURE [bar] 0\n    PRESSURE_TOLERANCE [bar] 100\n"
                 f"    OPTIMIZER LBFGS\n    MAX_ITER {MAX_ITER}\n    MAX_FORCE [eV*angstrom^-1] {MAX_FORCE}\n  &END CELL_OPT\n"
                 "  &PRINT\n    &RESTART\n      BACKUP_COPIES 3\n      &EACH\n        CELL_OPT 1\n      &END EACH\n"
                 "    &END RESTART\n  &END PRINT\n&END MOTION\n")
    return "\n".join(L)


# ----------------------------------------------------------------- output parsing
RE_E = re.compile(r"ENERGY\| Total FORCE_EVAL \( QS \) energy\s*\[(?:hartree|a\.u\.)\]:?\s+([-\d.Ee+]+)")
RE_GRID = re.compile(r"count for grid\s+(\d+):\s+(\d+)\s+cutoff")
RE_VOL = re.compile(r"CELL\| Volume \[angstrom\^3\]:\s+([-\d.Ee+]+)")
RE_NP = re.compile(r"Total number of message passing processes\s+(\d+)")
RE_PRESS = re.compile(r"OPT\| Pressure deviation \[bar\]\s+([-\d.Ee+]+)")


def last_segment(path):
    t = Path(path).read_text(errors="replace")
    i = t.rfind("PROGRAM STARTED AT")
    return t[i:] if i >= 0 else t


def is_done(d, proj, run_type):
    out = Path(d) / f"{proj}.out"
    if not out.exists():
        return False
    t = last_segment(out)
    if "PROGRAM ENDED AT" not in t or not RE_E.findall(t):
        return False
    if run_type == "ENERGY":      # grid scan: MAX_SCF 1 is never converged by design
        return True
    return "SCF run NOT converged" not in t and "OPTIMIZATION COMPLETED" in t


def parse_point(d, proj):
    t = last_segment(Path(d) / f"{proj}.out")
    e = RE_E.findall(t)
    vol = RE_VOL.findall(t)
    pr = RE_PRESS.findall(t)
    return dict(E_ha=float(e[-1]), volume=float(vol[-1]) if vol else None,
                pressure_dev_bar=float(pr[-1]) if pr else None,
                grids=[int(c) for _, c in RE_GRID.findall(t)[:5]],
                ranks=int(RE_NP.search(t).group(1)) if RE_NP.search(t) else None,
                opt_ok="OPTIMIZATION COMPLETED" in t)


# ----------------------------------------------------------------- job pool
class Pool:
    """Free-CPU pool: each job asks for n CPUs (default slot_cpus) and is pinned to them."""
    def __init__(self, ncpus, slot_cpus, threads=2):
        self.threads = threads
        cpus = sorted(os.sched_getaffinity(0))[:ncpus]
        if len(cpus) < ncpus:
            print(f"[WARN] only {len(cpus)} CPUs in affinity (asked {ncpus})", flush=True)
        self.slot_cpus = min(slot_cpus, len(cpus))
        self.free = cpus
        self.total = len(cpus)
        self.lock = threading.Condition()

    def run(self, d, proj, run_type, label, n=0):
        inp = Path(d) / f"{proj}.inp"
        rst = Path(d) / f"{proj}-1.restart"
        if is_done(d, proj, run_type):
            print(f"[skip] {label} already done", flush=True)
            return
        with self.lock:
            n = min(n or self.slot_cpus, self.slot_cpus_total())
            while len(self.free) < n:
                self.lock.wait()
            cpus, self.free = self.free[:n], self.free[n:]
        use = rst.name if rst.exists() and run_type == "CELL_OPT" else inp.name
        ranks = max(1, len(cpus) // self.threads)   # hybrid: harpia gives 900 MB/CPU, 1 GiB/rank measured
        cmd = ["taskset", "-c", ",".join(map(str, cpus)), "mpirun", "-n", str(ranks),
               "cp2k.psmp", "-i", use, "-o", f"{proj}.out"]
        env = dict(os.environ, OMP_NUM_THREADS=str(self.threads), MKL_NUM_THREADS=str(self.threads), OMP_PLACES="cores", I_MPI_PIN="0",
                   I_MPI_HYDRA_BOOTSTRAP="fork", FI_PROVIDER="shm", OMP_STACKSIZE="512M")
        t0 = time.time()
        print(f"[start] {label} on CPUs {cpus[0]}-{cpus[-1]} ({ranks} ranks x {self.threads} thr) input={use}", flush=True)
        with open(Path(d) / f"{proj}.stdout", "ab") as so:
            rc = subprocess.call(cmd, cwd=d, env=env, stdout=so, stderr=subprocess.STDOUT)
        ok = is_done(d, proj, run_type)
        print(f"[{'ok' if ok else 'FAIL'}] {label} rc={rc} {time.time() - t0:.0f}s", flush=True)
        with self.lock:
            self.free = sorted(self.free + cpus)
            self.lock.notify_all()

    def slot_cpus_total(self):
        return self.total

    def run_all(self, jobs):  # jobs: list of (dir, proj, run_type, label); concurrent
        ths = [threading.Thread(target=self.run, args=j) for j in jobs]
        for t in ths:
            t.start()
            time.sleep(1)
        for t in ths:
            t.join()


# ----------------------------------------------------------------- stages
def scan_grid(a, pool, struct, phase, values, cutoff, rel):
    lat, sp, frac = struct
    root = Path(a.workdir) / ("01_cutoff" if phase == "cutoff" else "02_relcutoff")
    jobs = []
    for v in values:
        c, r = (v, rel) if phase == "cutoff" else (cutoff, v)
        tag = f"{phase}{v:g}"
        d = root / tag
        d.mkdir(parents=True, exist_ok=True)
        proj = f"TiNbTa_{tag}"
        (d / f"{proj}.inp").write_text(make_input(proj, "ENERGY", lat, sp, frac, c, r))
        jobs.append((str(d), proj, "ENERGY", tag))
    if a.gen_only:
        return None
    pool.run_all(jobs[::-1])  # most expensive first
    rows = []
    for d, proj, _, tag in jobs:
        if not is_done(d, proj, "ENERGY"):
            print(f"[WARN] {tag}: not finished, excluded")
            continue
        p = parse_point(d, proj)
        p["value"] = float(tag.replace(phase, ""))
        rows.append(p)
    nat = len(sp)
    ref = rows[-1]["E_ha"]
    for r in rows:
        r["dE_mev_atom"] = abs(r["E_ha"] - ref) * HA_EV * 1000 / nat
        g = r["grids"]
        r["grid_ok"] = bool(g) and g[0] > 0 and sum(g[1:]) > g[0]
    rec = rec_energy_only = None
    for i, r in enumerate(rows):
        if all(x["dE_mev_atom"] < a.threshold for x in rows[i:]):
            rec_energy_only = rec_energy_only if rec_energy_only is not None else r["value"]
            if r["grid_ok"]:
                rec = r["value"]
                break
    note = ""
    if rec is None and rec_energy_only is not None:
        rec, note = rec_energy_only, "grid criterion (finest grid used, majority on coarse) NOT met; energy criterion only"
    print(f"\n{phase:>10} {'E [Ha]':>18} {'dE meV/at':>10}  Gaussians per grid (1=finest)")
    for r in rows:
        print(f"{r['value']:>10g} {r['E_ha']:>18.10f} {r['dE_mev_atom']:>10.3f}  {r['grids']} {'grid-ok' if r['grid_ok'] else ''}")
    print(f"=> recommended {phase.upper()} = {rec} Ry {note}\n", flush=True)
    (root / "summary.json").write_text(json.dumps(dict(phase=phase, recommended=rec, note=note, rows=rows,
                                                      threshold_mev_atom=a.threshold), indent=1))
    return rec


def stage_kmesh(a, pool, struct, cutoff, rel):
    lat, sp, frac = struct
    root = Path(a.workdir) / "03_kmesh"
    jobs = []
    for m in a.meshes:
        tag = "k%d-%d-%d" % m
        d = root / tag
        d.mkdir(parents=True, exist_ok=True)
        proj = f"TiNbTa_{tag}"
        (d / f"{proj}.inp").write_text(make_input(proj, "CELL_OPT", lat, sp, frac, cutoff, rel, kmesh=m, cell_ref=1.15))
        jobs.append([str(d), proj, "CELL_OPT", tag, 0])
    # CPUs per mesh ~ sqrt(#k-points), summing to the whole allocation
    w = [math.sqrt(m[0] * m[1] * m[2]) for m in a.meshes]
    n = [max(2, 2 * int(round(pool.total / 2 * x / sum(w)))) for x in w]   # even: 2 threads/rank
    while sum(n) > pool.total:
        n[n.index(max(n))] -= 2
    while sum(n) < pool.total:
        n[n.index(max(n))] += 2
    for j, k in zip(jobs, n):
        j[4] = k
    if a.gen_only:
        return None
    pool.run_all(jobs[::-1])
    rows = []
    for (d, proj, _, tag, _n), m in zip(jobs, a.meshes):
        if not is_done(d, proj, "CELL_OPT"):
            print(f"[WARN] {tag}: CELL_OPT not completed, excluded")
            continue
        p = parse_point(d, proj)
        p["mesh"] = list(m)
        rows.append(p)
    nat = len(sp)
    ref = rows[-1]
    for r in rows:
        r["dE_mev_atom"] = abs(r["E_ha"] - ref["E_ha"]) * HA_EV * 1000 / nat
        r["dV_pct"] = 100 * (r["volume"] - ref["volume"]) / ref["volume"]
    rec = None
    for i, r in enumerate(rows):
        if all(x["dE_mev_atom"] < a.threshold and abs(x["dV_pct"]) < a.vol_threshold for x in rows[i:]):
            rec = r["mesh"]
            break
    print(f"{'mesh':>8} {'E [Ha]':>18} {'dE meV/at':>10} {'V [A^3]':>10} {'dV %':>7} {'|P-P0| bar':>9}")
    for r in rows:
        print(f"{'x'.join(map(str, r['mesh'])):>8} {r['E_ha']:>18.10f} {r['dE_mev_atom']:>10.3f} "
              f"{r['volume']:>10.3f} {r['dV_pct']:>7.3f} {str(r['pressure_dev_bar']):>9}")
    print(f"=> recommended k-mesh = {rec} (dE < {a.threshold} meV/atom and |dV| < {a.vol_threshold} %)\n", flush=True)
    (root / "summary.json").write_text(json.dumps(dict(recommended=rec, rows=rows), indent=1))
    return tuple(rec) if rec else None


def parse_restart(path):
    t = Path(path).read_text()
    sub = t[t.upper().find("&SUBSYS"):]
    cell = re.search(r"&CELL\b(.*?)&END CELL\b(?!_)", sub, re.S).group(1)
    cell = re.sub(r"&CELL_REF.*?&END CELL_REF", "", cell, flags=re.S)
    lat = []
    for k in "ABC":
        m = re.search(rf"^\s*{k}\s+(?:\[[^\]]*\]\s+)?([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)", cell, re.M)
        lat.append([float(x) for x in m.groups()])
    co = re.search(r"&COORD\b(.*?)&END COORD", sub, re.S).group(1)
    scaled = bool(re.search(r"^\s*SCALED\s*(\.TRUE\.|T|TRUE)?\s*$", co, re.M | re.I))
    sp, xyz = [], []
    for ln in co.splitlines():
        p = ln.split()
        if len(p) >= 4 and p[0].upper() not in ("SCALED", "UNIT"):
            try:
                v = [float(x) for x in p[1:4]]
            except ValueError:
                continue
            sp.append(re.match(r"[A-Z][a-z]?", p[0]).group(0))
            xyz.append(v)
    if not scaled:   # Cartesian -> fractional
        inv = inv3(lat)
        xyz = [[sum(c[i] * inv[i][j] for i in range(3)) for j in range(3)] for c in xyz]
    return lat, sp, xyz


def inv3(m):
    d = det3(m)
    c = [[(m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3]
           - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]) / d for j in range(3)] for i in range(3)]
    return c


def stage_cellopt(a, pool, struct, cutoff, rel, mesh):
    lat, sp, frac = struct
    root = Path(a.workdir) / "04_cellopt"
    root.mkdir(parents=True, exist_ok=True)
    proj = "TiNbTa_bulk"
    (root / f"{proj}.inp").write_text(make_input(proj, "CELL_OPT", lat, sp, frac, cutoff, rel, kmesh=mesh, cell_ref=1.15))
    if a.gen_only:
        return
    pool.run_all([(str(root), proj, "CELL_OPT", "bulk CELL_OPT")])
    if not is_done(root, proj, "CELL_OPT"):
        sys.exit("[FAIL] production CELL_OPT did not complete (see 04_cellopt/TiNbTa_bulk.out)")
    p = parse_point(root, proj)
    rl, rsp, rfr = parse_restart(root / f"{proj}-1.restart")
    write_poscar(root / "bulk_relaxed.vasp", rl, rsp, rfr, "Nb4Ta1Ti15 CP2K PBE relaxed")
    v0 = abs(det3(struct[0]))
    res = dict(natoms=len(sp), E_ha=p["E_ha"], E_ev=p["E_ha"] * HA_EV, E_ev_per_atom=p["E_ha"] * HA_EV / len(sp),
               volume_A3=abs(det3(rl)), volume_vasp_A3=v0, dV_vs_vasp_pct=100 * (abs(det3(rl)) - v0) / v0,
               pressure_dev_bar=p["pressure_dev_bar"], cutoff=cutoff, rel_cutoff=rel, kmesh=list(mesh),
               lattice=rl, abc=[norm(v) for v in rl], species=rsp, frac=rfr,
               a_bcc_equiv_A=(2 * abs(det3(rl)) / len(sp)) ** (1 / 3))
    (root / "bulk.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in res if k not in ("lattice", "species", "frac")}, indent=1), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--structure", default=str(HERE / "VASP_references" / "CONTCAR (1)"))
    p.add_argument("--workdir", default=str(HERE / "bulk"))
    p.add_argument("--stage", default="all", choices=["all", "cutoff", "relcutoff", "kmesh", "cellopt"])
    p.add_argument("--ncpus", type=int, default=48)
    p.add_argument("--threads", type=int, default=2, help="OpenMP threads per MPI rank")
    p.add_argument("--slot-cpus", type=int, default=0, help="ranks per concurrent point (default: per stage)")
    p.add_argument("--cutoff-values", type=float, nargs="+", default=[300, 400, 500, 600, 700, 800, 900, 1000])
    p.add_argument("--rel-values", type=float, nargs="+", default=[20, 30, 40, 50, 60, 70, 80, 100])
    p.add_argument("--rel-start", type=float, default=60, help="REL_CUTOFF during the CUTOFF scan")
    p.add_argument("--cutoff", type=float, help="force CUTOFF (skip/override stage result)")
    p.add_argument("--rel-cutoff", type=float, help="force REL_CUTOFF")
    p.add_argument("--kmesh", type=int, nargs=3, help="force k-mesh")
    p.add_argument("--mesh-list", default=";".join(",".join(map(str, m)) for m in MESHES))
    p.add_argument("--threshold", type=float, default=1.0, help="meV/atom")
    p.add_argument("--vol-threshold", type=float, default=0.5, help="|dV| %% for k-mesh")
    p.add_argument("--gen-only", action="store_true", help="write inputs only")
    a = p.parse_args()
    a.meshes = [tuple(int(x) for x in t.split(",")) for t in a.mesh_list.split(";")]
    Path(a.workdir).mkdir(parents=True, exist_ok=True)
    struct = read_poscar(a.structure)
    print(f"structure: {len(struct[1])} atoms, V={abs(det3(struct[0])):.3f} A^3, "
          f"valence e-={sum(GTH_Q[s] for s in struct[1])}", flush=True)
    dec_path = Path(a.workdir) / "decisions.json"
    dec = json.loads(dec_path.read_text()) if dec_path.exists() else {}
    run = lambda s: a.stage in ("all", s)

    def pool(slot):
        return Pool(a.ncpus, a.slot_cpus or slot, a.threads)

    cutoff, rel = a.cutoff or dec.get("cutoff"), a.rel_cutoff or dec.get("rel_cutoff")
    mesh = tuple(a.kmesh) if a.kmesh else (tuple(dec["kmesh"]) if dec.get("kmesh") else None)
    if run("cutoff") and not a.cutoff:
        cutoff = scan_grid(a, pool(6), struct, "cutoff", a.cutoff_values, None, a.rel_start)
        dec["cutoff"] = cutoff
    if run("relcutoff") and not a.rel_cutoff:
        if cutoff is None and not a.gen_only:
            sys.exit("no CUTOFF available: run stage cutoff or pass --cutoff")
        rel = scan_grid(a, pool(6), struct, "relcutoff", a.rel_values, cutoff or max(a.cutoff_values), a.rel_start)
        dec["rel_cutoff"] = rel
    dec_path.write_text(json.dumps(dec, indent=1))
    if run("kmesh") and not a.kmesh:
        if a.gen_only:
            cutoff, rel = cutoff or 600, rel or 50
        elif cutoff is None or rel is None:
            sys.exit("need CUTOFF and REL_CUTOFF before kmesh")
        mesh = stage_kmesh(a, pool(a.ncpus), struct, cutoff, rel)
        dec["kmesh"] = list(mesh) if mesh else None
    dec_path.write_text(json.dumps(dec, indent=1))
    if run("cellopt"):
        if a.gen_only:
            cutoff, rel, mesh = cutoff or 600, rel or 50, mesh or (5, 4, 4)
        elif None in (cutoff, rel, mesh):
            sys.exit(f"incomplete decisions: cutoff={cutoff} rel={rel} mesh={mesh} -- extend the scans")
        stage_cellopt(a, pool(a.ncpus), struct, cutoff, rel, mesh)
    print("DONE", dec, flush=True)


if __name__ == "__main__":
    main()
