# SECOND_BRAIN: contexto denso do projeto CP2K (`~/work_cp2k`)

> Arquivo de trabalho para as IAs. O cofre (`~/cofre/02_Memory/cp2k/*.md`) guarda o resumo e
> aponta para cá. **Reconstruído em 2026-09-24 na máquina `access`.** O `SECOND_BRAIN.md` original
> ficou em `coaraci:~/work_cp2k` (sem versionamento, inacessível nesta data). As seções 6 e
> 9b–9e do original (RI-HFXk) estão resumidas aqui a partir da nota do cofre `cp2k/tio2.md`.

## 0. Máquinas e sincronização
- `access` (GridUNESP, login `access`): 56 nós de 56 cores, 126 GB por nó; partições `short`
  (1 d), `medium` (7 d), `long` (30 d); cerca de 730 jobs pendentes por partição em 2026-09-24.
  Módulos: `cp2k/2023.1`(default), `cp2k/2026.1`, `cp2k/2026.1.gpu`.
- `coaraci` (UNICAMP): onde o pipeline foi reescrito pela primeira vez (2026-09-21/22, CP2K
  2026.2). SSH deu timeout em 2026-09-24, e essa versão nunca foi versionada.
- Repo: `git@github.com:robertoramosjr/cp2k_scripts.git`. **Só código e docs** (`.gitignore` em
  whitelist). Dados (`ZrO2/`, `TiO2_smoke/`) nunca entram.
- Python: `~/.conda/envs/cp2k_env` (pymatgen 2026.7.31). O `vasp_env` tem pymatgen 2025.

## 1. Fatos da build `cp2k/2026.1` (access), verificados
- Linkada contra o **MPICH 4.3.2 do toolchain**. Com `srun` puro, sobem N cópias de 1 rank.
  Usar `mpiexec -n $SLURM_NTASKS -bind-to core:$OMP_NUM_THREADS` (pinning testado: rank0→0,1;
  rank1→2,3).
- O módulo exporta `CP2K_DATA`, mas o binário lê `CP2K_DATA_DIR`: sempre
  `export CP2K_DATA_DIR="$CP2K_DATA"`.
- **SCF não convergido aborta o run** (`IGNORE_CONVERGENCE_FAILURE` é False por padrão). Só o
  scan de grid (MAX_SCF 1) liga essa keyword.
- `&OT` + qualquer `&KPOINTS` (inclusive 1 1 1) aborta. Em Γ, omitir a seção.
- PDOS: só `&DFT%PRINT%PDOS` (COMPONENTS, NLUMO, APPEND). `&DFT%PRINT%DOS` só tem DOS total
  (DELTA_E); `&DOS%PDOS` é rejeitado.
- `&HF%RI` tem as keywords RI_FLAVOR, RI_METRIC, EPS_FILTER, MEMORY_CUT; **sem**
  `KP_RI_MEMORY_ESTIMATE` (essa é da 2026.2). A sintaxe RI-HFX + k-points passa no `--check`.
- Não existe a subseção `&XC_FUNCTIONAL/&LIBXC`: cada funcional libxc é uma subseção própria
  (ex. `&HYB_GGA_XC_HSE06`).
- Schema completo: `cp2k.psmp --xml` gera `cp2k_input.xml` (34 MB), que pode ser consultado com
  ElementTree (ver `scripts/probe_build.sh`).
- Bases: DZVP-MOLOPT-SR-GTH-q12 existe para Zr e Ti (26 funções esféricas); O q6 tem 13. TZVP e
  TZV2P-MOLOPT-SR de Zr estão em `BASIS_MOLOPT_UCL`.

## 2. Pipeline (ver CP2K_PIPELINE.md)
`01_grid_convergence → 02_kmesh_convergence → 03_cellopt → {04_bands, 05_pdos} → 06_slab_cut → 07_slab_opt`,
com parsers em `parsers/`, jobs em `jobs/`, `core/cp2k_blocks.py` como fonte única e `--basis`
obrigatório.
Validação real em 2026-09-24 (login, 2 ranks, ZrO₂ cúbico, cutoff 200–300): `run_dir` do
`common.sh` rodou com 2 ranks reais, e o 01 mais o `parse_grid` funcionaram; ver seção 5 para
bandas/PDOS.

## 3. Pegadinhas (as antigas continuam válidas)
- `PROGRAM ENDED AT` não significa convergência. Os jobs checam `OPTIMIZATION COMPLETED`.
- `sort` em arrays tem que ser `sort -V`. Não editar um script de array com tarefas pendentes.
- O writer binário de restart de k-points (`.kp`) crasha sob I/O concorrente. Nos arrays,
  `&SCF%PRINT%RESTART OFF`.
- pymatgen: `Slab.is_polar()` é sempre False sem estados de oxidação. O 06 atribui
  `--oxidation Zr:4,O:-2`.
- pymatgen `get_tasker2_slabs()` falha no t-(110); o 06 tem uma reconstrução Tasker III→II
  própria (`tasker3_to_2`).
- Com argparse, `--planes -1,1,1` é lido como opção. Usar `--planes="-1,1,1;..."`.

## 4. RI-HFXk (resumo do coaraci, TiO2 rutilo): pausado
Na 2026.2 funcionou, mas a memória é inviável (38–77 GiB/rank em EPS_PGF_ORB 1e-3/1e-4 para 6
átomos em 2×2×2). Tudo foi testado e fechado com runs reais: `RI_FLAVOR MO` aborta (é
incompatível com `&SMEAR`); `EPS_STORAGE_SCALING` não muda nada; `EPS_FILTER` dá só ~20%; Γ-only
com R_c menor custa `erfc(ω·R_c)` da troca (30% perdidos em 2,9 Å); diagonalização sem RI dá
segfault ("Only RI-HFX is implemented for K-points"). O usuário foi para o **CRYSTAL** para
bandas híbridas. Não re-testar.

## 5. Projetos
### ZrO2: slabs PBE → gamma → Wulff → mapa de morfologia (ATIVO)
- Pasta: `ZrO2/PAW_PBE/` (mapa em `ZrO2/PAW_PBE/README.md`). VASP original em `PAW_PBE/legacy/`.
- Bulks de partida: `CONTCAR` do Etot VASP (ISIF=2 sobre a rede convergida por Stress): c a=5.115;
  t a=3.621, c=5.278; m a=5.186, b=5.243, c=5.374, β=99.63° (P2₁/c padrão).
- Planos (literatura, ver o README do ZrO2): m (-111), (111), (-101), (011), (110); t (101), (111),
  (001), (100), (110); c (111), (110), (100), e **2 a definir com o usuário**.
- Status 2026-09-24: 01 (cutoff) gerado para as 3 fases e slabs cortados. **Nada submetido.**

### TiO2 (`TiO2_smoke/`, pipeline antigo): referência histórica. Os tempos valem para 1 rank × 2 threads.
