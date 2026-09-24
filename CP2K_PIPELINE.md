# CP2K_PIPELINE: diagnóstico e desenho do pipeline CP2K

> **Proveniência.** Este documento foi **reconstruído em 2026-09-24 na máquina `access`**. O
> original (`CP2K_PIPELINE.md` / `cp2k_pipeline.md`) e a primeira implementação
> (`core/cp2k_blocks.py`, `00_basis_probe` … `06_bands_hybrid`, `SECOND_BRAIN.md`) foram escritos
> em 2026-09-21/22 em `coaraci:~/work_cp2k`. Nada disso foi versionado: `work_cp2k` não era repo
> git, e `coaraci` estava inacessível (timeout de SSH) nesta data. A reconstrução parte de três
> fontes: o prompt de reescrita, as notas do cofre (`02_Memory/cp2k/tio2.md`) e os problemas
> reais encontrados no pipeline antigo desta máquina (agora em `legacy/`). **Quando o `coaraci`
> voltar, comparar os dois lados e mesclar** (ver PARTE 3, item 6). O SECOND_BRAIN do coaraci
> chegou em 2026-09-24 ([SECOND_BRAIN_COARACI.md](SECOND_BRAIN_COARACI.md)); as divergências
> estão em SECOND_BRAIN.md §5. A spec original ainda não chegou.

---

## PARTE 1: diagnóstico do pipeline antigo (`legacy/scripts`, `legacy/cp2k_fluxo.md`)

1.1 **MPI.** O `cp2k/2026.1` usa o MPICH 4.3.2 do próprio toolchain. Os jobs carregavam
    `intel/mpi/2017` e chamavam `srun` sem plugin PMI, o que subia N cópias de 1 rank do mesmo
    cálculo, todas escrevendo no mesmo `.out` (`Total number of message passing processes 1`).
    Todos os tempos do `TiO2_smoke` são de 1 rank × 2 threads.
1.2 **`&KPOINTS` em Γ.** O gerador sempre escrevia `&KPOINTS` (até `1 1 1`), o que coloca o CP2K
    no caminho de k-points, onde `&OT` aborta ("OT not possible with kpoint calculations").
1.3 **Drift entre scripts.** `cutoff_convergence.py` e `kmesh_stress_convergence.py` passavam
    `mixing_alpha` para um `build_scf_block` que não aceitava o argumento (TypeError). Cada
    gerador duplicava o template do input.
1.4 **`ADDED_MOS` derivado do número de elétrons** (cópia do `NBANDS=NELECT` do VASP): era caro e
    não tinha motivo físico no CP2K. O certo é um valor fixo pequeno, default 10.
1.5 **Etapas redundantes.** Um `CELL_OPT` (k-mesh) era seguido de `GEO_OPT` de produção e de um
    "verify-relax". O `CELL_OPT` já relaxa íons e célula juntos: basta **um** passo de produção.
1.6 **Convergência de grid 1D.** Varria só `CUTOFF` com SCF completo. O método oficial do CP2K
    (how-to "Converging the CUTOFF and REL_CUTOFF") é 2D:
    - `RUN_TYPE ENERGY`, Γ-only, **`MAX_SCF 1`** a partir do mesmo chute `ATOMIC`. O erro de SCF
      é idêntico em todos os pontos e cancela nas diferenças, que então medem só o erro de grid.
      Na build 2026.1 o SCF não convergido **aborta** por padrão, então é obrigatório
      `IGNORE_CONVERGENCE_FAILURE TRUE` (só neste passo).
    - `PRINT_LEVEL MEDIUM` imprime `count for grid N`: quantas Gaussianas foram para cada nível do
      multigrid.
    - Fase 1: varre `CUTOFF` com `REL_CUTOFF` fixo (60). Fase 2: varre `REL_CUTOFF` no `CUTOFF`
      escolhido.
    - Critério: o menor valor com ΔE < limiar (meV/átomo) contra o maior valor testado, para ele
      e para todos os maiores. Além disso, o grid mais fino precisa ser usado (contagem > 0),
      mas com a **maioria das Gaussianas nos grids grossos**.
1.7 **Smearing ausente em k-points:** as SCFs de óxidos com k-mesh denso oscilavam. Agora `&SMEAR`
    FERMI_DIRAC 300 K entra sempre que houver k-points, com `ADDED_MOS` para acomodá-lo.
1.8 **PDOS/bandas no mesmo run, com k-mesh da relaxação multiplicado.** Agora o PDOS usa uma
    supercélula em Γ com OT (barato e sem ambiguidade), e o broadening é feito em Python.
    As bandas usam `kpath.prim` do pymatgen, com a célula do input igual à do caminho.
1.9 `&SCF%PRINT%RESTART` implícito. Agora é sempre explícito: ON na produção, OFF em arrays
    concorrentes (o writer binário `.kp` crashou sob I/O concorrente, com SIGABRT em
    `write_kpoints_restart`).

## PARTE 2: desenho novo

```
core/cp2k_blocks.py    blocos de input (fonte única) + DEFAULTS + CLI comum (--basis obrigatório)
core/cp2k_output.py    leitores de output (energia, grids, ranks, convergência)
01_grid_convergence.py --phase cutoff|relcutoff   (1.6)
02_kmesh_convergence.py CELL_OPT completo por malha; critério duplo ΔE meV/át + ΔV %
03_cellopt.py          CELL_OPT de produção, passo único: CELL_REF 1.15x, STRESS_TENSOR ANALYTICAL,
                       BFGS/LBFGS, EXTERNAL_PRESSURE 0, MAX_FORCE em eV/Å, &SCF%PRINT%RESTART ON
04_bands.py            ENERGY + &PRINT%BAND_STRUCTURE, k-mesh do 03, ADDED_MOS 10, kpath.prim
05_pdos.py             supercélula (--supercell NX NY NZ), Γ + OT, &PRINT%PDOS COMPONENTS, NLUMO -1
06_slab_cut.py         slabs simétricos, estequiométricos e apolares (Tasker III→II), ordenados
                       por ligações quebradas por área
07_slab_opt.py         GEO_OPT do slab reescalado na rede do 03, c ∥ z, correção de dipolo
parsers/parse_{grid,kmesh,cellopt,bands,pdos,surface_energy}.py
jobs/common.sh, job_scan_array.sh (01/02/lote de slabs), job_single.sh (03/04/05/slab)
scripts/probe_build.sh o que a build aceita (PARTE 3, itens 1 e 2)
```

Regras do `core/cp2k_blocks.py`:
- Γ-only → `&OT` (DIIS, FULL_SINGLE_INVERSE) + `&OUTER_SCF`, **sem** a seção `&KPOINTS`.
- k ≠ Γ → `&DIAGONALIZATION` (STANDARD, `EPS_ADAPT 0.01`) + `&MIXING` (Broyden, α 0.4) +
  `&SMEAR ON` (FERMI_DIRAC, 300 K) + `ADDED_MOS 10` (nunca derivado do número de elétrons).
- `&QS EXTRAPOLATION USE_GUESS`, `EPS_DEFAULT 1E-12`; `&KPOINTS PARALLEL_GROUP_SIZE -1`.
- CELL_OPT: `EXTERNAL_PRESSURE [bar] 0` explícito (o default do CP2K é 100 bar). Otimizador
  BFGS até 10 átomos, LBFGS acima.
- `build_cell_block(structure, cell_ref_factor)` escreve `&CELL_REF` (fator linear nos vetores)
  só em `CELL_OPT`.
- Tudo o que é tunável vem de `DEFAULTS` via argparse. Os únicos valores fixos ficam em `jobs/*.sh`.

Jobs: `-N 1 --ntasks-per-node=14 --cpus-per-task=2 --mem-per-cpu=2194M`, `module cp2k/2026.1`,
`CP2K_DATA_DIR="$CP2K_DATA"`, lançamento com **`mpiexec -n $SLURM_NTASKS -bind-to core:2`** (MPICH
Hydra; nunca `srun` puro, ver 1.1) e aviso quando o número de ranks do CP2K ≠ `SLURM_NTASKS`.
São restart-safe via `<project>-1.restart`. Sucesso significa `OPTIMIZATION COMPLETED` em
GEO/CELL_OPT, e `PROGRAM ENDED` com energia em ENERGY. Nunca basta ter terminado.

Convenções: argparse em tudo, scripts em inglês, rodam de qualquer diretório e agem no diretório
corrente ou no `--output-dir`, cada etapa grava `run_meta.json`/`scan_meta.json` para a próxima.

## PARTE 3: pendências

1. **Interface DOS/PDOS.** `probe_build.sh` em 2026-09-24 com `cp2k/2026.1`: `&DFT%PRINT%PDOS`
   existe (COMPONENTS, NLUMO); `&DFT%PRINT%DOS` existe só com DOS total (DELTA_E), **sem
   `&PDOS` aninhado** (o `--check` rejeita). Default `--dos-interface legacy`. O layout
   unificado era o da 2026.2 no coaraci.
2. **RI-HFX com k-points.** A sintaxe (`&HF%RI` + `&KPOINTS`) passa no `--check` da 2026.1. O
   teste real (`probe_build.sh --run`, só via sbatch) ainda não foi feito aqui. No coaraci, com
   a 2026.2, funcionou, mas a memória ficou inviável: `KP_RI_MEMORY_ESTIMATE` deu 38–77 GiB/rank
   para o TiO2 de 6 átomos. Essa keyword **não existe** na 2026.1.
3. `parse_bands.py` foi escrito contra o formato documentado do `.bs`. Validar no primeiro run
   real (feito com ZrO₂ cúbico em 2026-09-24, ver SECOND_BRAIN.md).
4. Bandas/DOS híbridas (`06_bands_hybrid` do coaraci) não foram portadas: o usuário foi para o
   CRYSTAL nesse ponto.
5. Espessura dos slabs: 06 corta com `--min-slab 12` Å. Falta testar a convergência de gamma com
   a espessura (+1 e +2 camadas nos planos de menor gamma).
6. Reconciliar com `coaraci:~/work_cp2k` quando o SSH voltar (diff de `core/`, `0?_*.py`, docs).
