# SECOND_BRAIN.md: unified reference for Claude (read FIRST after any context compaction)

Audience: Claude, future sessions. Dense; update §10 (STATE) at the end of every work block.
**Unified 2026-09-24 on `access`** from the coaraci SECOND_BRAIN (2026-09-21/22, originally in English, kept
VERBATIM below where it is the primary source: rules, coaraci machine, verified 2025.1/2026.2 facts, HSE06/RI-HFXk,
TiO2 state) and the access notes (2026-09-24). Originals: `legacy/docs/SECOND_BRAIN_COARACI.md`,
`legacy/docs/SECOND_BRAIN_ACCESS_2026-09-24.md`. Spec: [CP2K_PIPELINE.md](CP2K_PIPELINE.md) (also unified).
Long coaraci session notes (only on coaraci): `~/.claude/projects/-home-users-roarj-work-cp2k/memory/*.md` (index MEMORY.md).
Cross-project brain: `~/cofre` (git; `git pull` at start, push at end; notes `02_Memory/cp2k/{tio2,zro2}.md`).

---------------------------------------------------------------------------------------------------
## 0. Rules of engagement
From coaraci (verbatim; machine-specific items marked where access differs):
- Chat with the user in **Portuguese**; scripts, comments, argparse help, docstrings in **English**.
- `cp2k --check` proves SYNTAX only. Never say "works" without a real run (use queue `teste` here).
- **MaxRSS at OOM time is a snapshot when killed = lower bound**, not the requirement.
- CP2K **appends** to an existing `-o` file: parsers/monitors must look only at the last `PROGRAM STARTED AT` segment
  (old ABORTs from earlier runs caused false alarms twice).
- On coaraci `srun` **without `--mpi=pmix`** runs each rank as an independent singleton
  (`GLOBAL| Total number of message passing processes 1`). Always verify that line == requested tasks.
- `/tmp` is **node-local** on compute nodes: sbatch in/out files must live on the shared FS (`~/work_cp2k/...`).
  The scratchpad (`/tmp/claude-1010076/.../scratchpad`) is login-node only (fine for XML dumps, helper scripts).
- Bash tool: no `sleep N` chains (blocked). Wait with `until <cond>; do sleep 5; done` (run_in_background) or Monitor.
  Monitor scripts must dedupe (only new lines) or they re-emit forever. Stop monitors when the user says not to poll.
- **The Write tool's PreToolUse hook has timed out twice (nothing written)**: fall back to `cat > file << 'EOF'` via Bash and verify with `ls`.
- Primary working dir shifts between tool calls: **use absolute paths**.
- Long unattended drivers: launch with `setsid nohup ... > log 2>&1 < /dev/null & disown`. A plain `nohup ... &` child got SIGINT (exit -2) mid-run once.
- Ask before big/expensive/irreversible actions; say plainly when something can't be done as asked
  (e.g. queue `short` and `cp2k/2026.1` do not exist on this machine).
- Don't spend the shared queue on runs I predict will fail: cheap feasibility probe first.
- The user is loose with names ("fila short", "CP2K_PIPELINE.md" = `cp2k_pipeline.md`). Check reality, say so briefly.
- The user accepted being corrected on premises ("wannier ran fine" was false; "small system => cheap" is false for RI-HFXk).

Added on access:
- On **access** the launcher is MPICH `mpiexec -n $SLURM_NTASKS -bind-to core:$OMP_NUM_THREADS` (NOT `srun --mpi=pmix`;
  cp2k/2026.1 is linked to its own toolchain MPICH 4.3.2). Queues `short/medium/long` DO exist here; `teste` does not.
- The user's rule (cofre): report cost BEFORE submitting; nothing heavy on the login node (a 2-rank, few-minute test is ok).
- pymatgen gotchas for surfaces: `Slab.is_polar()` needs oxidation states; `SlabGenerator` and `Slab(...)` rotate the lattice
  unless `reorient_lattice=False`; `--planes -1,1,1` is eaten by argparse (use `--planes="-1,1,1;..."`).
- Never publish secrets: the legacy scripts had a hardcoded Materials Project API key (removed before the first push).

---------------------------------------------------------------------------------------------------
## 1. User
Roberto (roarj), computational materials scientist (Unicamp). VASP background (ISIF, NBANDS, ICHARG=11 mental model).
CP2K for oxides; TiO2 now; production targets 80-120 atoms (e.g. BaZrS3). VS Code remote. Values: doc-verified claims, evidence,
being told when a premise is off. Sometimes grants autonomy ("faça como achar melhor"); otherwise wants to choose on trade-offs.
- On access the account is `rramos` (GridUNESP); on coaraci `roarj`. Also runs VASP (`~/work_vasp`) and Turbomole projects.

---------------------------------------------------------------------------------------------------
## 2. Machines
### access (GridUNESP login `access`) — PRIMARY since 2026-09-24; this is the "target cluster" the coaraci notes talk about
- 56 nodes x 56 cores, 126000 MB/node; partitions `short` (1 d, default), `medium` (7 d), `long` (30 d), `gpu`; ~730 jobs pending per partition (2026-09-24).
- Modules: `cp2k/2023.1` (default), `cp2k/2026.1` (psmp+ssmp, GCC 14.3 + MPICH 4.3.2 toolchain), `cp2k/2026.1.gpu`. The module sets `CP2K_DATA`,
  the binary reads `CP2K_DATA_DIR` -> `export CP2K_DATA_DIR="$CP2K_DATA"`.
- Python: `~/.conda/envs/cp2k_env` (pymatgen 2026.7.31, matplotlib); `vasp_env` has pymatgen 2025. SSH to coaraci times out from here.
- Repo: `~/work_cp2k` = `git@github.com:robertoramosjr/cp2k_scripts.git` (`main`), whitelist `.gitignore` (data never versioned).
### coaraci ( `coaraci.ifi.unicamp.br`; Slurm 23.11.4, Oracle Linux 8)
- Nodes: 48 cores (2x24, no SMT), 120400 MB. Partitions: `paralela`, `par48-x/i` (single node), `par480-x/i` (multi-node, MinNodes=2,
  per-user job limits), `serial`, `gpu-*`, `clauber`, **`teste`** (2 nodes, 30 min — my test bed), **`fat`** (1 node, 2 TB, 3 d; AllowAccounts
  restricted; was 100% allocated with a queue on 2026-09-21 -> not used). **No `short/medium/long`.**
- CP2K modules: `cp2k/2024.1-...`, `cp2k/2025.1-gcc-12.2.0-fod7ldk` (only one with `cp2k.psmp`; openmpi 5.0.3 + pmix 5.0.1),
  `cp2k/2025.2-gcc-14.2.0-sjpfo5a` (default, BROKEN: aborts even on `--xml`). No 2026.x anywhere (latest public release is 2025.2).
- `CP2K_DATA` NOT set by the module. Data dir: `/opt/spack/opt/spack/linux-oracle8-zen2/gcc-12.2.0/cp2k-2025.1/share/cp2k/data`.
- Non-interactive shells: `source /etc/profile.d/*lmod* ; module load ...`.
- Slurm >= 22.05: export `SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK` else ranks get 1 core.
- Default `python3` has pymatgen, numpy 2.x (`np.trapz` removed -> `np.trapezoid`), matplotlib.
- Old smoke data: `~/work_cp2k/TiO2_smoke/` (24-atom TiO2; `kmesh_stress_convergence/` holds REAL converged CELL_OPT outputs = good parser test data).
### Local CP2K 2026.2 (installed 2026-09-21 by me, from the user's upload `cp2k-2026.2-Linux-gnu-x86_64.ssmp`)
- `~/software/cp2k-2026.2/bin/cp2k.ssmp` (428 MB, STATIC, prebuilt; rev 67b5da876d, GCC 15.3; flags omp libint fftw3 libxc libxsmm spglib dftd4 tblite ...). **ssmp = OpenMP only, NO MPI, one node**
  (=> 1 process x N threads; KP_NGROUPS must stay 1; the singleton-MPI trap does not apply). `source ~/software/cp2k-2026.2/env.sh` sets PATH, CP2K_DATA_DIR, OMP_STACKSIZE.
- Data dir = symlink to the cluster's cp2k/2025.1 data (the binary embeds a nonexistent build path). Verified: same input gives the same energy as 2025.1 to 4e-13 Ha.
- Run: `sbatch -p par48-x|teste -N1 -n1 -c48 --mem=115G ... --export=ALL,OMP_NUM_THREADS=48 --wrap='source ~/software/cp2k-2026.2/env.sh; cp2k.ssmp -i X.inp -o X.out'`.
  With the pipeline jobs: `--export=ALL,CP2K_MODULE=,CP2K_EXE=$HOME/software/cp2k-2026.2/bin/cp2k.ssmp,CP2K_DATA=...` (jobs now skip `module load` when CP2K_MODULE is empty), `-n1 -c48`.
- 2026.2 schema (checked): has `&DFT%PRINT%DOS` with `CURVE/PDOS/LDOS/R_LDOS` (spec pendency 1 = YES); `SORT_BASIS EXP`; `ADMM_TYPE` (ADMMS ok); `EPS_DIIS`; `WFN_RESTART_FILE_NAME` (in &DFT);
  `&HF/&RI/&PRINT` logical keys `KP_RI_MEMORY_ESTIMATE` and `KP_RI_PROGRESS_BAR` ("rough UPPER BOUND of memory for an RI-HFXk ENERGY calc; needs a fair amount of computing first").
  None of these exist in 2025.1 (SORT_BASIS, RI print keys) -> only validate them with the 2026.2 binary.

---------------------------------------------------------------------------------------------------
## 3. Repo layout (`~/work_cp2k` on access = cp2k_scripts)
```
CP2K_PIPELINE.md (+ symlink cp2k_pipeline.md, what the cofre's _Sources/cp2k points to)   SECOND_BRAIN.md   README.md
core/cp2k_blocks.py     all input blocks + DEFAULTS + add_common_args (--basis mandatory) + structure/.restart readers
core/cp2k_output.py     output readers (last PROGRAM STARTED segment only): energy, grid counts, ranks, convergence
01_grid_convergence.py  --phase cutoff|relcutoff (2D method)      02_kmesh_convergence.py  CELL_OPT per mesh (--k-list|--densities)
03_cellopt.py           production CELL_OPT (CELL_REF 1.15x, EXTERNAL_PRESSURE 0, BFGS<=10 at/LBFGS, --keep-space-group)
04_bands.py             ENERGY + BAND_STRUCTURE on kpath.prim     05_pdos.py  Gamma supercell OT PDOS (--dos-interface legacy|unified)
06_slab_cut.py          symmetric/stoichiometric/non-polar slabs (Tasker III->II), ranked by broken bonds, [uvw] + hkl self-check
07_slab_opt.py          slab on the CP2K lattice, c||z, GEO_OPT or ENERGY (--vacuum, --run-type)
parsers/parse_{grid,kmesh,cellopt,bands,pdos,surface_energy,vacuum}.py
jobs/common.sh (env, run_cp2k via mpiexec, is_done, run_dir with restart) + job_scan_array.sh + job_single.sh
scripts/probe_build.sh  schema + --check probes (PDOS interface, RI-HFX+k); --run only inside sbatch
legacy/scripts, legacy/cp2k_fluxo.md   old access pipeline;   legacy/docs/   original coaraci + access docs
ZrO2/ TiO2_smoke/       data (not versioned)
```
**Only on coaraci (not ported):** `00_basis_probe.py`, `06_bands_hybrid.py` (-> port as `08_bands_hybrid.py`), BASES registry
(`--basis dzvp-molopt-sr|ccgrb-d|ccgrb-t|pob-tzvp-rev2`, GAPW, ADMM aux), `scripts/run_smoke.py`, parsers `--json`,
`check_restart_matches_input`, `[max_iter]` detection, `&HF/&MEMORY MAX_MEMORY`, `--functional PBE0|HSE06`, `smoke/` data.
Access `--basis` takes the raw family name (e.g. `DZVP-MOLOPT-SR-GTH`) and appends `-qN`.

---------------------------------------------------------------------------------------------------
## 4. Spec digest: see CP2K_PIPELINE.md (PART 1 = findings 1.1-1.19, PART 2 = pipeline, PART 6 = surfaces)
Divergences between the two implementations, deliberately kept (flag if the user objects): access `NGRIDS 5` (coaraci 4);
access `&OUTER_SCF` on the OT branch (coaraci none); access `MAX_SCF 300` (coaraci 100 PBE / 300 hybrid; the rutile PDOS needed 150);
grid criterion in meV/atom (coaraci 1e-5 Ha total). Adopted from coaraci on 2026-09-24: `EXTERNAL_PRESSURE 0`, LBFGS > 10 atoms,
last-segment parsing, unified-DOS layout (NLUMO at `&DOS`), 2026.2 `.pdos` header parsing, `KEEP_SPACE_GROUP` option.

---------------------------------------------------------------------------------------------------
## 5. Verified CP2K facts
### 5a. cp2k/2026.1 on access (2026-09-24)
- Default aborts on unconverged SCF (`IGNORE_CONVERGENCE_FAILURE` F). `&OT` + any `&KPOINTS` aborts. `EXTERNAL_PRESSURE` default 100 bar.
- DOS: `&DFT%PRINT%PDOS` (COMPONENTS, NLUMO, APPEND) exists; `&DFT%PRINT%DOS` = total DOS only (DELTA_E); `&DOS%PDOS` rejected by --check.
- `&HF%RI`: RI_FLAVOR, RI_METRIC, EPS_FILTER, MEMORY_CUT; **no KP_RI_MEMORY_ESTIMATE**. RI-HFX + k-points passes --check (no real run yet).
- `MOTION/CELL_OPT` and `GEO_OPT` have `KEEP_SPACE_GROUP`, `KEEP_SYMMETRY` (CELL_OPT), `EPS_SYMMETRY`. No `&XC_FUNCTIONAL/&LIBXC` wrapper.
- Bases: Zr/Ti DZVP-MOLOPT-SR-GTH-q12 = 26 spherical functions, O q6 = 13; Zr TZVP/TZV2P-MOLOPT-SR in BASIS_MOLOPT_UCL.
- Real 2-rank login runs (c-ZrO2): mpiexec gives 2 ranks; `.bs` format = coaraci's; PBE c-ZrO2 indirect gap X->Gamma 3.30 eV (200-300 Ry);
  2026.1 `.pdos` header = one line (kind + E(Fermi)); c-ZrO2 at REL_CUTOFF 60 / 200-300 Ry: most Gaussians still on grid 1.
- Cost calibration (access TiO2_smoke, 24 atoms, 72 k, 350 Ry): ~30 s per SCF step on 2 cores (k-point diag dominated); legacy runs = 1 rank x 2 threads.
### 5b. From coaraci (2025.1 unless noted; verbatim)
#### k-point restart (verified 2026-09-21 with 2026.2 ssmp, PBE rutile 2x2x2)
- `&SCF%PRINT%RESTART ON` writes `<project>-RESTART.kp` for k-points (`.wfn` is Gamma-only); read back with `SCF_GUESS RESTART` + `&DFT WFN_RESTART_FILE_NAME x-RESTART.kp`
  (output shows `read_kpoints_restart`). A .kp written with `&QS EPS_PGF_ORB 1e-12` was read by a run with the default (1e-6): 2 SCF steps instead of 15, same energy/bands.
- 2026.2 schema (`cp2k.ssmp --xml`, paths WITHOUT the CP2K_INPUT root, e.g. `/FORCE_EVAL/DFT/XC/HF/RI`): RI keys `KP_NGROUPS`(1), `KP_STACK_SIZE`(16 in 2026.2, 32 in 2025.1), `EPS_PGF_ORB`(1e-5),
  `MEMORY_CUT`(3); `&RI/&PRINT` LOGICAL keys `KP_RI_MEMORY_ESTIMATE`, `KP_RI_PROGRESS_BAR`; `&DFT` `SORT_BASIS EXP`, `AUTO_BASIS RI_HFX MEDIUM`; `ADMM_TYPE ADMMS` (sets METHOD/purification/scaling
  itself: do not also write METHOD/ADMM_PURIFICATION_METHOD); `&SCF EPS_DIIS`. Estimate keyword "is more accurate when restarting from a (GGA) wavefunction; forces need more memory".
- RI-HFXk truncation: "a sphere of radius R_c must fit in the BvK supercell" (cell x k-mesh). rutile 4x4x4 with R_c 6 A: c-axis 4x2.96=11.8 A < 12 -> `check_truncation_radius` warns (needs 5 points along c).
- RI-HFXk prints `KP-HFX_RI_INFO|` lines (parallel groups, stack size, RI extension radius, unit cells covered, avg sgf in extended RI bases, image cells considered): what parse_basis_probe reads.
#### Input schema (query `cp2k.psmp --xml` with ElementTree; XML lands in cwd as cp2k_input.xml)
- `CELL_REF` is a subsection of `&CELL`. CELL_OPT `EXTERNAL_PRESSURE` default 100 bar; `PRESSURE_TOLERANCE` 100 bar.
- `MAX_SCF 1` aborts unless `IGNORE_CONVERGENCE_FAILURE T`. `&SMEAR ON`, `EXTRAPOLATION USE_GUESS` valid.
- `KPOINT_SET NPOINTS n` -> n+1 points (endpoints inclusive); one set per path segment works; `UNITS B_VECTOR` default.
- `&DFT%PRINT%DOS` exists (only `&EACH`); `&CURVE` does NOT (new DOS interface unavailable). `KP_RI_MEMORY_ESTIMATE` does NOT exist.
- `KP_NGROUPS` docstring: N subgroups => ~Nx speed AT THE COST OF Nx MEMORY. `MEMORY_CUT` 3->5 gave nothing.
#### Output formats
- Energy `ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]  x` (older `[a.u.]`). Atoms: `- Atoms:  N`.
- Grid: under `MULTIGRID INFO`: `count for grid N: COUNT cutoff [a.u.] X` (grid 1 = FINEST) + `total gridlevel count`.
- CELL_OPT per step: `OPT| Step number`, `OPT| Total energy [hartree]`, `OPT| Internal pressure [bar]`, `CELL| Volume [angstrom^3]`; done marker
  `GEOMETRY OPTIMIZATION COMPLETED` (substring `OPTIMIZATION COMPLETED`), then "Reevaluating energy at the minimum" with the final CELL block+energy.
  `<proj>-1.restart` (+ `-1_N.restart`, `.bak-N`) written each step.
- Restart nests `&CELL_REF` INSIDE `&CELL` (naive "last A/B/C" returns the 1.15x cell). Restart keeps `MAX_ITER` + `STEP_START_VAL`: a resume after MAX_ITER
  does nothing until MAX_ITER is raised IN THE RESTART (jobs print `[max_iter]`; verified the sed fix then continues).
- `.bs`: `# Set n: 2 special points, k k-points, B bands` / `#  Special point i kx ky kz NAME` / `#  Point n Spin s: kx ky kz w` / `#   Band Energy [eV] Occupation` rows.
  `.pdos`: header `# Projected DOS for atomic kind Ti ... E(Fermi) = x a.u.`, cols `MO Eigenvalue[a.u.] Occupation s py pz px d-2..d+2 (f..)`, files `<proj>-k1-1.pdos`.
- SCF: `SCF run converged in N steps` / `SCF run NOT converged`; end `PROGRAM ENDED AT`; Fermi: `Fermi energy:` (Ha).
#### Physics sanity (rutile 6 atoms)
- PBE 2x2x2 cutoff 400: direct gap at Gamma ~1.9 eV (lit ~1.8); PDOS VB=O2p, CB=Ti3d. Grid scan at REL_CUTOFF 60: finest grid still ~60% of Gaussians at
  400 Ry -> parse_grid correctly refuses to recommend on a short scan.

---------------------------------------------------------------------------------------------------
## 6. HSE06 / RI-HFXk (coaraci, verbatim) — PAUSED; do not re-derive, do not re-propose closed options
#### Required input (implemented in core, `--functional HSE06`)
- XC `&HYB_GGA_XC_HSE06` + `&HF` FRACTION 0.25, `&INTERACTION_POTENTIAL` SHORTRANGE OMEGA 0.11, `&SCREENING` EPS_SCHWARZ 1e-6 SCREEN_ON_INITIAL_P TRUE.
- Basis files BASIS_MOLOPT + `BASIS_ADMM` + `BASIS_ADMM_MOLOPT`; ADMM `AUX_FIT`: **O `cFIT3`** (BASIS_ADMM), **Ti `cFIT10`** (BASIS_ADMM_MOLOPT)
  (names are library names, not `cFIT3-PBE-qN` — that guess was legacy bug #1). `RI_AUX` optional (absence accepted, no memory effect).
  Generators warn if CP2K_DATA_DIR/CP2K_DATA data files lack the aux name; override with `--kind-overrides {"Ti":{"aux_basis_set":"..."}}`.
- **k != Gamma needs RI-HFX**: `&HF/&RI` with `EPS_PGF_ORB = sqrt(EPS_DEFAULT)` (default RI 1e-5 != QS 1e-6 aborts) and ADMM purification **NONE**
  (MO_DIAG needs OT=Gamma; NONE_DM rejected). Gamma: plain HFX + ADMM MO_DIAG with OT. (4 legacy input bugs, each visible only after fixing the previous.)
#### Cost — the real problem (never assume "small system => cheap")
- RI-HFXk cost is set by the OMEGA/EPS_SCHWARZ radius (~14 A), not atom count. Extended RI basis: TiO2 24 atoms 39 cells/470 image cells; **Si2 2x2x2: 2193 image
  cells, >8 GB, still running at 900 s (4 cores, 120 GB)**; **rutile 6 atoms 2x2x2: 1294 image cells, 5709 sgf, OOM in ~1 min on a 112 GB node (12 ranks x 4 cpus,
  ~9.4 GB/rank at kill) — both for bands (ENERGY) and CELL_OPT (so HFX stress support is untested)**. 2x2x2 is already the smallest sensible mesh.
- Legacy 24-atom TiO2: 1 node 12x4 ~10-12 GB/rank OOM at 5-12 min right after `KP-HFX_RI_INFO`; multi-node made per-rank memory ~2x WORSE; MEMORY_CUT, RI_AUX tier
  down/removed, k-mesh 18x8x8->14x6x6 (KP-HFX_RI_INFO identical: OMEGA-driven), 24x2 split: nothing closed it. Last legacy submission (par480-x, 3 nodes, job 1786122)
  was left queued; outcome unknown. Untested lever: loosen EPS_SCHWARZ (real accuracy tradeoff -> needs user OK). Rejected: KP_NGROUPS (more memory).
- Gamma HSE06 (rutile 6 atoms, OT): works, 89 s for 100 SCF steps, but the OT gradient was only 5e-6 at step 100 -> hybrids get MAX_SCF 300. CP2K WARNS
  "cutoff radius larger than half the minimal cell dimension -> may lead to unphysical total energies" (small cells + SHORTRANGE, Gamma): PDOS-with-HSE06 on small
  supercells is physically dubious; use larger supercells for real work.
- Literature: VASP HSE06 on TiO2 uses coarse k (2x2x2..7x7x7). CP2K's official RI-HFXk example (graphene) uses TRUNCATED R=5 A + KP_NGROUPS 16 + ADMM.
- Practical route (standard): relax with PBE, hybrid only for the electronic structure where feasible.

#### 2026-09-21 investigation "why can't we do HSE06 bands" (user's premise "HSE06 = GGA + plane waves, should be cheap" is wrong: HSE06 is a HYBRID; in CP2K orbitals are
Gaussians (GPW) and exact exchange needs 4-centre integrals with no reciprocal-space FFT shortcut; VASP does HFX in k-space with plane waves. Cost model differs fundamentally.)
- Docs (manual ri_kpoints + Bussy&Hutter JCP 160, 064116 (2024)): RI-HFXk is built around a TRUNCATED Coulomb operator, R_c=5 A ("converged from ~6 A"); "a limited range potential is required";
  "a sphere of radius R_c must fit in the BvK supercell" (k-mesh*cell > 2*R_c per direction). Designed for "small unit cells and dense k-meshes". EPS_PGF_ORB default 1e-5 "accurate and fast"
  (the pipeline used 1e-6 via EPS_DEFAULT 1e-12 -> costlier; RI EPS_PGF_ORB must equal sqrt(QS EPS_DEFAULT): use `--eps-default 1e-10` for hybrids).
- Facts measured on rutile 6 atoms, k=2x2x2, 12 ranks x 4 cpus, 112 GB, CP2K 2025.1 (all OOM within 1-2 min, BEFORE the first SCF iteration):
  HSE06 SHORTRANGE (auto radius 14.28 A): 1294 image cells (eps 1e-12) / 1168 (eps 1e-10). PBE0 TRUNCATED R_c=5: 295 images; R_c=6: 345 images. All die at ~10-12 GB/rank.
- **`POTENTIAL_TYPE MIX_CL_TRUNC` (HSE06 = 1/r - erf(wr)/r truncated) ABORTS with k-points: "HFX potential not available for K-points (NYI)"**. For k-points only SHORTRANGE (auto radius) and TRUNCATED exist.
- Knobs that did NOT rescue it: R_c 5/6 A, `KP_STACK_SIZE` 1 and 4 (`--kp-stack-size`), `AUTO_BASIS RI_HFX SMALL` (no effect: with ADMM the RI basis is derived from the AUX_FIT basis,
  printed "RI HFX Basis Set Ti-RI-AUX-cFIT10 = 133 fns, O-RI-AUX-cFIT3 = 49 fns"; cFIT10/cFIT3 already smallest), 4 ranks x 12 threads (31.5 GB/rank, same total), minimal SZV orbital basis
  (identical image cells / "4253 sgf in extended RI bases", same OOM => cost sits in the RI/image-cell part, NOT the orbital basis). `RI_AUX` declarations are irrelevant to HFX (HFX uses type RI_HFX).
- True memory requirement is UNKNOWN (only lower bounds >115 GB). Ways to learn it: `&HF/&RI/&PRINT KP_RI_MEMORY_ESTIMATE` (documented for trunk; NOT in 2025.1 schema; probably in 2026.1 —
  try it there; syntax per docs: `&PRINT` containing the bare key `KP_RI_MEMORY_ESTIMATE`), or a bigger node (`fat`, 2 TB: contended, needs user OK).
- Feasible today in CP2K: Gamma-only hybrid on supercells (standard HFX+ADMM, works, 12 atoms in ~1 min) — but needs cell >= 2*R_c (HSE06 14 A => impossible for small cells; PBE0-TC R_c 4-6 A => supercell >= ~9-12 A).
  Hybrid *gap at Gamma* is obtainable that way (rutile gap is direct at Gamma); full hybrid bands are not, on <=120 GB nodes.
- Honest framing for the user (wants to migrate from VASP): for small-cell hybrid band structures plane-wave codes are far cheaper; CP2K's edge is large cells at Gamma. Practical CP2K route: PBE bands + hybrid Gamma-supercell
  gap/DOS (scissor-type correction), or RI-HFXk on a big-memory/multi-node allocation after getting the memory estimate.
- Pipeline options added for this: `--functional PBE0`, `--hfx-potential {SHORTRANGE,TRUNCATED,MIX_CL_TRUNC}`, `--hfx-cutoff-radius`, `--kp-stack-size`, `--auto-basis-ri-hfx`. Experiment dir `smoke/hyb/`
  (inputs hA..hD, s1..s3, t1; slurm outputs; jobids*.txt).

#### 2026-09-21 (later) 2026.2 measurements — the memory estimator WORKS and EPS_PGF_ORB is THE lever (`smoke/eps_sweep/`, `smoke/probe_build_run/out26*`)
- Local CP2K 2026.2 = `~/software/cp2k-2026.2/bin/cp2k.ssmp` (ONE rank, OpenMP only => KP_NGROUPS must stay 1; source `env.sh`). `KP-HFX_RI_INFO| Estimated peak memory usage per MPI rank (MiB)` is printed
  right after the image-cell info, BEFORE the first SCF step, provided the setup itself still fits in memory (`&HF/&RI/&PRINT KP_RI_MEMORY_ESTIMATE T`). Sample per-rank estimates (1 rank, 48 threads):
  rutile 6 at, DZVP, 2x2x2, R_c 2.9 A: EPS_PGF_ORB 1e-3 -> 38.2 GiB (ext RI 1794 sgf, 141 image cells; 1st SCF step 535 s); 1e-4 -> 77.2 GiB (2508 sgf, 212 cells; still no 2nd step in 10 min).
  1e-5 (4253 sgf, 251 cells), 1e-6 (5709 sgf, 299 cells) and the real probe 4x4x4/R_c 6/1e-6 (427 cells) DIE (OOM) before printing => > ~115 GiB.
  Memory ~ (extended-RI sgf)^2 (38 -> 77 GiB for sgf 1794 -> 2508: x2.02 vs x1.95 predicted) => extrapolated ~215 GiB at 1e-5 and ~390 GiB at 1e-6 for this trivial 2x2x2 case (EXTRAPOLATION, not measured).
  Si2 (2 atoms) 2x2x2, eps 1e-6: HSE06 SHORTRANGE (no cutoff, 24.4 A, 2193 image cells) 42.5 GiB vs PBE0 TRUNCATED R_c 3.5 (446 cells) 11.8 GiB.
- The R_c effect on image cells is real but the extension radius (~5.7 A) comes from EPS_PGF_ORB + basis; ccGRB-D needs 12856 sgf in the extended RI basis vs 5709 for DZVP-MOLOPT-SR at 1e-6 (2x).
- PBE-level effect of EPS_PGF_ORB (relaxed rutile 6x6x6, `smoke/eps_pbe/`): 1e-12/1e-6/1e-5/1e-4 all give gap 1.7706 eV, dE(1e-4 vs 1e-12) = 9e-7 Ha; 1e-3 shifts gap by 1.4 meV, dE 3.4e-5 Ha.
  Whether RI-HFXk tolerates 1e-4 (RI default is 1e-5) is UNVERIFIED for the hybrid gap. The pipeline default stays sqrt(EPS_DEFAULT)=1e-6 (`--eps-pgf-orb` to change; QS/RI kept identical automatically).
- `--hfx-potential` now only SHORTRANGE|TRUNCATED (MIX_CL_TRUNC removed, it aborts with k-points). PBE0 default R_c 6.0 A (was 5.0).

#### coaraci 9c. 2026-09-22 (later) — untried RI-HFXk manual alternatives, actually tested now (`smoke/ri_flavor_test/`)
User asked whether manual-documented HFX memory alternatives had been tried before concluding CP2K can't do it here. Full `&HF/&RI` schema (2026.2) dumped and cross-checked against everything tried before (SECOND_BRAIN sect. 6): found 3 real untried keywords, tested all 3 on the same baseline (rutile 6 atoms, 2x2x2, R_c 2.9 A, EPS_PGF_ORB 1e-3, dzvp-molopt-sr, RHO default = 39324 MiB/rank estimate):
- **`RI_FLAVOR MO`** (vs default RHO; docs say "RHO typically requires more memory"): **ABORTS** — "RI_FLAVOR MO is not consistent with smearing. Please use RI_FLAVOR RHO" (hfx_types.F:1144). Our pipeline always uses `&SMEAR` on the k-point branch (cp2k_pipeline.md 1.3, required for the Fermi energy) → **hard-blocked, not usable here, not just "didn't help."**
- **`EPS_STORAGE_SCALING`** (0.01 → 1.0, storage threshold for 3-center integrals): 39324 → 39044 MiB, **negligible (~0.7%)**.
- **`EPS_FILTER`** (1e-9 → 1e-6, DBT tensor contraction filter): 39324 → 31424 MiB, **real, ~20% reduction** — the one genuine, previously-untried lever that helps. NOT free: loosens the tensor filtering precision, an actual accuracy tradeoff (untested how it propagates to the SCF/gap), and still an order of magnitude short of what production EPS_PGF_ORB (1e-5/1e-6, extrapolated 215-390 GiB) would need even with a 20% cut applied (still ~170-310 GiB vs ~120 GiB/node available).
**Conclusion given to the user**: the manual's alternatives were NOT all tried before this — asking was the right call. Now they are (RI_FLAVOR MO ruled out structurally, EPS_STORAGE_SCALING negligible, EPS_FILTER real-but-insufficient-alone). Nothing found today closes the order-of-magnitude gap; `EPS_FILTER` loosening is worth stacking with everything else (EPS_PGF_ORB, KP_NGROUPS=1, R_c) if a bigger-memory node ever becomes available, but does not make it feasible on a 120 GB node by itself. Inputs/outputs in `smoke/ri_flavor_test/` (`flavor_rho`, `flavor_mo`, `eps_filter_loose`, `eps_storage_loose`).

#### coaraci 9d. 2026-09-22 (later still) — "traditional DIAGONALIZATION for hybrid k-point bands, no RI" empirically re-confirmed impossible (`smoke/diag_no_ri_test/`)
User proposed skipping RI-HFX entirely: hybrid functional + real k-mesh + plain `&SCF/&DIAGONALIZATION` (SCF_GUESS ATOMIC, ADDED_MOS, no OT) with NO `&HF/&RI` section -- i.e. the snippet's DIAGONALIZATION/ADDED_MOS/SCF_GUESS ATOMIC part is exactly what 04_bands.py/06_bands_hybrid.py already do for any k-point run (OT cannot do k-points, not new), but the snippet had no `&HF/&RI` at all. Built and RAN exactly that (Si2, HSE06, 2x2x2 k-mesh, &HF without &RI) on 2026.2: got past atomic-guess setup, entered "SCF WAVEFUNCTION OPTIMIZATION", then **segfaults with zero diagnostic** the instant the first HFX build starts (no clean [ABORT] box this time, unlike the 2025.1-era "Only RI-HFX is implemented for K-points" message documented in [[cp2k-hybrid-kpoints-ri-fix]] hfx_admm_utils.F -- same underlying restriction, cruder failure mode on this build). **Conclusion given to the user: this is not a real alternative** -- DIAGONALIZATION is already mandatory and already used (not something new that helps), and CP2K has NO non-RI implementation of hybrid exchange under k-point sampling at all; removing &HF/&RI doesn't reduce memory, it just removes the only code path capable of doing the calculation, so it crashes before doing any useful work. Do not re-propose "skip RI, use plain HFX with k-points" without new CP2K version evidence.

#### coaraci 9e. 2026-09-22 (final) — USER PAUSED CP2K, moving to CRYSTAL
After exhausting every real avenue for hybrid-functional band structures on this hardware (sect. 9b-9d: RI-HFXk memory wall ~38->390 GiB depending on EPS_PGF_ORB, RI_FLAVOR MO blocked by smearing, EPS_FILTER only -20%, non-RI diagonalization segfaults, Gamma-only truncated-R_c route quantified as a real ~13-30% exchange error unless the supercell grows into the hundreds of atoms), the user decided: **CP2K work is paused here for now**, moving to CRYSTAL (LCAO code with a more mature periodic exact-exchange implementation, better suited to hybrid+k-point bands than CP2K's newer RI-HFXk) for this specific need (hybrid band structures); VASP remains their tool for what VASP already does well. Not abandoned permanently -- if resumed, read sect. 6 and 9b-9d FIRST (all the "did you try X" alternatives already have real, tested answers; don't re-derive).
**State to resume from**: PBE full pipeline (00-05) works end to end, real numbers in sect. 9b. Hybrid bands (06_bands_hybrid.py) implemented and passes cp2k.ssmp --check + a minimal feature test, but never converged on the real TiO2 system -- memory, not a code bug. GGA+U (DFT_PLUS_U, per-KIND, no RI/HFX dependency, should NOT hit the same wall) was proposed to the user as an untested-but-promising alternative for the TiO2 gap; user declined for now (going to CRYSTAL instead of pursuing it in CP2K).
**Untouched**: no CP2K jobs left running (checked before pausing); `~/cofre` (the multi-project vault, see sect. 3) has this project's `02_Memory/cp2k/tio2.md` updated to reflect the pause.

---------------------------------------------------------------------------------------------------
## 7. Job scripts
### access (jobs/, current)
- `common.sh`: module cp2k/2026.1, CP2K_DATA_DIR, OMP_*; `run_cp2k` = mpiexec + rank check; `is_done` by RUN_TYPE (OPT banner / ENDED+energy);
  `run_dir` resumes from `<project>-1.restart` (renames the old .out), writes `status.txt` (OK / IN_PROGRESS / FAILED). Does NOT yet detect `[max_iter]`.
- `job_scan_array.sh DIR` (one task per subdir, `sort -V`), `job_single.sh DIR`. Header: short, 1 d, `-N1 --ntasks-per-node=14 --cpus-per-task=2 --mem-per-cpu=2194M`.
- Under SLURM `$0` is a spool copy: `common.sh` is sourced from `~/work_cp2k/jobs/` as fallback.
### coaraci (verbatim, for porting)
- Environment block overridable by env: `CP2K_MODULE` (default cp2k/2026.1), `CP2K_EXE`, `SRUN_MPI_FLAG` (default empty), `OPT_DONE_MARKER`; exports CP2K_DATA_DIR
  from CP2K_DATA if set, OMP_*, SRUN_CPUS_PER_TASK, `ulimit -s unlimited`; cd to submit dir. No `--time` set (target limits unknown).
- Every job runs `check_mpi` (warns if CP2K's MPI process count != SLURM_NTASKS). grid: loops `*.inp`, skips finished. kmesh: loops `k*/` (array index supported).
  cellopt: `run_cellopt` -> `[skip]/[start]/[resume]/[done]/[incomplete](exit 1)/[max_iter](exit 3)/[FAILED](exit 2)`. bands/pdos: single ENERGY run, skip if done, check `.bs`/`*.pdos`.
- **Run on coaraci**: `sbatch -p teste -N1 --ntasks-per-node=12 --cpus-per-task=4 --mem-per-cpu=2400M --time=00:29:00
  --export=ALL,CP2K_MODULE=cp2k/2025.1-gcc-12.2.0-fod7ldk,CP2K_DATA=<data dir>,SRUN_MPI_FLAG=--mpi=pmix jobs/job_X.sh <arg>` (`--mem-per-cpu` overrides the header; `--mem=0` conflicts).
- `CP2K_MODULE` uses `${CP2K_MODULE-default}` (NO colon): set-but-empty means "no module load, use CP2K_EXE as is" (with `:-` an empty value silently fell back to the default; fixed 2026-09-21).
- job_bands_hybrid.sh (also runs the 00 probes): validates `SLURM_NTASKS % KP_NGROUPS` (KP_NGROUPS read from the input, default 1) -> `[ABORT]` exit 3; requires a `.kp` guess when the input has WFN_RESTART_FILE_NAME.
- The jobs were generated from one template (scratchpad `make_jobs.py`, may be gone): edit the .sh files directly and keep the shared functions identical.

---------------------------------------------------------------------------------------------------
## 8. Projects
### ZrO2 slabs (ACTIVE, access) — `~/work_cp2k/ZrO2/PAW_PBE/README.md` is the map; cofre `cp2k/zro2.md`
- VASP originals in `PAW_PBE/legacy/` (untouched). Starting bulks = VASP Etot CONTCARs (c a=5.115; t a=3.621 c=5.278; m 5.186/5.243/5.374, beta 99.63, P2_1/c std).
- Planes (literature): m (-111),(111),(-101),(011),(110); t (101),(111),(001),(100),(110); c (111),(110),(100) + 2 TBD with the user.
- 27 slabs cut (06), all verified independently (plane via [uvw], symmetric, stoichiometric, zero dipole, no duplicates); `_1` = fewest broken bonds.
- Legacy VASP slabs rejected (wire c-100_1, non-stoichiometric, asymmetric/polar, duplicates, missing m-(-111)/(-101)); VASP mono slab_1 never converged SCF (tilted normal + IDIPOL=3).
- Cost model: 27 slabs ~39 node-h (8-117), `_1` only ~16 (3-49); x1.6 with 28-core tasks; all fit `short`.
### TiO2 rutile (coaraci) — PAUSED 2026-09-22 (user went to CRYSTAL for hybrid bands). Last states, verbatim:
#### TiO2 STATE (coaraci) — 2026-09-22 (HFX &MEMORY / PDOS-interface check)
- **User asked**: verify a PBE0 HFX snippet found online (`&HF/&MEMORY MAX_MEMORY`, `EPS_STORAGE_FILE`, `T_C_G_DATA`) against what the pipeline does, suspecting a missing piece. Also: add findings to `cp2k_pipeline.md` (the "manual" —
  clarified the user meant this file, not an external PDF). Result documented in `cp2k_pipeline.md` PARTE 5 and memory `cp2k_hfx_memory_and_pdos_interface.md`: the snippet's `EPS_STORAGE_FILE` keyword does not exist (CP2K aborts;
  real name `EPS_STORAGE`); `T_C_G_DATA` was never missing (resolves via CP2K_DATA_DIR like other data files); `&HF/&MEMORY MAX_MEMORY` (default 512 MiB/rank, stingy) WAS a real gap but only in the Gamma-only conventional-HFX
  path (`build_xc_block`, `05_pdos.py` with PBE0/HSE06 at Gamma) — fixed, now always written, default 4096, `--hfx-max-memory`. Confirmed this is UNRELATED to the RI-HFXk k-point memory wall (different CP2K subsection/code path).
- **Side discovery while regression-testing 05_pdos.py end-to-end on 2026.2**: `&DFT%PRINT%PDOS` standalone no longer exists there (moved inside a new `&DFT%PRINT%DOS` wrapper, NLUMO moved up) — exactly the cp2k_pipeline.md 1.9/
  PARTE-3-item-1 risk, now confirmed. Fixed: `build_pdos_block(dos_interface="legacy"|"unified")`, `05_pdos.py --dos-interface` (default legacy, unchanged behaviour/regression). Also fixed `parsers/parse_pdos.py`: the .pdos header
  itself changed on 2026.2 (kind/E(Fermi) split across lines, always full s/p/d/f columns) — `read_pdos()` now scans all header lines instead of assuming fixed positions; verified against a real 2026.2 file and a synthetic legacy one.
- **New, NOT yet fixed, found only while verifying the above**: a plain PBE Gamma-only OT PDOS run on rutile (`05_pdos.py`, no supercell, CUTOFF 900/REL_CUTOFF 40, fresh ATOMIC guess) did NOT converge in the default MAX_SCF 100
  (gradient 4.65e-6 at step 100, EPS_SCF default 1e-6) — needed `--max-scf 150`. Independent of dos-interface (reproduced with both). Never previously verified end-to-end on this exact combination; flagged, not fixed, out of scope.
- Regression re-verified after all of today's edits: 01/02/03/05(legacy PBE) generator outputs are still byte-identical to before this session; only Gamma-hybrid HSE06/PBE0 output changed (adds &MEMORY, intentional).

#### TiO2 prior STATE (coaraci) — 2026-09-21 (after the hybrid-bands addendum)
- **PBE full scan on par48-x: DONE** (`smoke/full_pbe/`, state `smoke_state.json`, logs `smoke/full_pbe_{b,c}.log`; driver finished, nothing running): CUTOFF 900 Ry / REL_CUTOFF 40 Ry (grid tolerance 1e-5 Ha = MY choice,
  the manual's 1e-8 Ha is unattainable: ~4e-6 Ha noise), k-mesh 6x6x6 (scan 2..10; picked with 1 meV/atom + 0.5 % V criteria, not at the scan edge), CELL_OPT 11 steps (4 min): a=b=4.641 A, c=2.969 A,
  V 62.449 -> 63.952 A^3 (+2.4 %), P 17 bar, E -180.9092784 Ha; PBE bands (04, 36 s, 12 ranks): direct gap at Gamma 1.7706 eV (VBM/CBM both at Gamma), `bands.bs/.csv/.png`, `bands-RESTART.kp` (67 MB, EPS_PGF_ORB 1e-12).
  Resume recipe if ever needed: `setsid nohup python3 scripts/run_smoke.py --cif smoke/rutile.cif --workdir smoke/full_pbe --basis dzvp-molopt-sr --start-at <stage> --stop-after bands --accept-fallback --grid-tol-energy 1e-5
  --sbatch-args "-p par48-x -N1 --ntasks-per-node=12 --cpus-per-task=4 --mem-per-cpu=2400M --time=08:00:00" --env CP2K_MODULE=cp2k/2025.1-gcc-12.2.0-fod7ldk --env SRUN_MPI_FLAG=--mpi=pmix --env CP2K_DATA=<2025.1 data dir> < /dev/null & disown`.
- **Addendum (bands with hybrid functional): items 1-7 IMPLEMENTED** (see sect. 3/5/7). Verified: PBE stage inputs byte-identical to pre-refactor (except 04: EPS_PGF_ORB + restart), `cp2k.ssmp --check` accepts 06/probe inputs,
  the .kp write/read cycle, job validation (ranks % KP_NGROUPS, missing/.wfn restart) with a fake srun, error paths (eps mismatch, 04 HSE06, .wfn, mismatched PBE mesh). NOT verified: any hybrid run to convergence
  (memory wall, sect. 6), GAPW+RI-HFXk beyond setup, ADMMS SCF behaviour, `.kp` guess reading in a hybrid run.
- **Numbers to give the user**: hybrid RI-HFXk memory is the blocker: measured estimates 38/77 GiB per rank at EPS_PGF_ORB 1e-3/1e-4 (rutile 2x2x2), > 115 GiB at 1e-5/1e-6; extrapolation ~215/~390 GiB. Suggest the user runs
  `00_basis_probe.py` + `KP_RI_MEMORY_ESTIMATE` on the 2026.x cluster and tries `--eps-pgf-orb 1e-4` or 1e-5 (validate the gap against the 1e-6 result if it ever fits).
- rutile 4x4x4 + PBE0 R_c 6 violates the BvK-fit rule along c (11.8 < 12 A): `check_truncation_radius` warns; use 4x4x5 or R_c <= 5.9.
- probe_build.sh on 2026.2 (`smoke/probe_build_run/out26*`, Si2 2x2x2): unified `&DFT%PRINT%DOS%CURVE` ACCEPTED (2026.2 has the new DOS interface). PBE0 RI-HFXk with the 06 keywords (ADMMS, SORT_BASIS, AUTO_BASIS,
  KP_RI_*, PROBE_EPS_PGF_ORB 1e-3, R_c 3.5) = PASS (normal termination; MAX_SCF 3 so NOT converged -- feature test only), est. 2.8 GiB/rank, 120-185 s per SCF step on 48 threads. HSE06 SHORTRANGE: ENGAGED only
  (1371 image cells at 1e-3, est. 12.4 GiB, no SCF step finished in 440 s): not proven to complete. New env vars: PROBE_EPS_PGF_ORB (default 1e-4), PROBE_RC, PROBE_VARIANTS, PROBE_TIMEOUT.
- Tell the user: on-disk `cp2k_pipeline.md` lacks sections 1.12 and "etapa 3c" (told once); line 201 still says `.wfn`; legacy `TiO2_smoke/dos_bands/job.sh` still passes a `.wfn` (old, untouched); `fat` earliest start 2026-10-20.
- Open questions for the user: EXTERNAL_PRESSURE default 0 OK?; ADDED_MOS 20 empty bands OK?; SRUN_MPI_FLAG on the new cluster?; move the 3 old scripts in `scripts/` to `legacy/`?; which `--alt-basis` to probe first?
- Untested: `--mp-id` path; 2026.1 cluster; 2026.1 output formats; HFX + analytical stress; Gamma hybrid with the new R_c 6.0 default on small supercells.

---------------------------------------------------------------------------------------------------
## 9. Timeline
- Aug 18-Sep 2: fixed legacy generators for HSE06+k-points (4 input bugs), found the srun singleton-MPI bug and SRUN_CPUS_PER_TASK gotcha, chased RI-HFXk memory (sect. 6).
  wannier_test (Aug 31) never reached Wannier90 (invalid KP_RI_MEMORY_ESTIMATE keyword, then OOM).
- Sep 21: user supplied `cp2k_pipeline.md`; rewrote everything, legacy -> `legacy/`, tested every stage with real runs (PBE rutile) on `teste`; bugs found by tests and fixed
  (kpoint scheme/parallel_group_size not propagated; probe verdict wording; `[max_iter]` detection via restart counters).
- Sep 21 (later): user asked (1) smoke re-run "on queue short with HSE06" and (2) this file. Added HSE06 to core, ran feasibility experiments (sect. 6), wrote `scripts/run_smoke.py`,
  `--json` in parsers, hybrid MAX_SCF default. `short` does not exist here -> used `teste` + CP2K 2025.1 as stand-in.
- Sep 22 (coaraci): hybrid alternatives closed (§6), CP2K paused.
- Sep 24 (access): coaraci unreachable; old access pipeline audited (srun/MPICH, KPOINTS-Gamma, drift bugs), pipeline rebuilt from the rewrite prompt,
  versioned as cp2k_scripts; ZrO2 PAW_PBE restructured, slabs audited and recut; user supplied the coaraci SECOND_BRAIN and spec -> both unified here.

---------------------------------------------------------------------------------------------------
## 10. STATE (update me) — 2026-09-24, access
- **Running/pending:** array `5209911` (ZrO2-vac, short, 65 tasks %20): vacuum convergence, single points (07 `--run-type ENERGY`) of the 13 `_1` slabs x
  vacuum 8/10/12/15/20 A, provisional 600/60 Ry, k-density 30, VASP lattice. Dir `ZrO2/PAW_PBE/vacuum_convergence/`. Parse:
  `python ~/work_cp2k/parsers/parse_vacuum.py --scan-dir ZrO2/PAW_PBE/vacuum_convergence` (threshold 5 mJ/m2). FIRST real SLURM use of jobs/common.sh:
  check `status.txt`, the `.err` for the rank warning, and that `message passing processes` = 28.
- **Ready, NOT submitted (needs user OK):** 01 cutoff scans (300-1000 Ry) in `ZrO2_*_pbe_cp2k/01_grid_convergence/cutoff/`.
- **Pending decisions (user):** 2 extra cubic planes; all 27 terminations or only `_1`; thickness convergence; cubic may relax to tetragonal in PBE.
- **Open code items:** `[max_iter]` detection in common.sh; ports from coaraci (§3); probe_build `--run` via sbatch.
