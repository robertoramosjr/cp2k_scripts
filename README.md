# cp2k_scripts

Pipeline CP2K (GPW, PBE) para sólidos e superfícies, usado no cluster GridUNESP. A especificação
e o diagnóstico estão em [CP2K_PIPELINE.md](CP2K_PIPELINE.md); o contexto denso (fatos da build,
pegadinhas, status dos projetos) em [SECOND_BRAIN.md](SECOND_BRAIN.md).

```
01_grid_convergence.py   CUTOFF / REL_CUTOFF (método 2D oficial)
02_kmesh_convergence.py  k-mesh com CELL_OPT completo (ΔE + ΔV)
03_cellopt.py            bulk de produção (CELL_OPT, passo único)
04_bands.py              bandas PBE (caminho do HighSymmKpath)
05_pdos.py               PDOS em supercélula Γ + OT
06_slab_cut.py           slabs simétricos, estequiométricos e apolares
07_slab_opt.py           GEO_OPT do slab na rede do CP2K
core/                    blocos de input + leitores de output
parsers/                 um parser por etapa (+ energia de superfície e Wulff)
jobs/                    SLURM (common.sh, job_scan_array.sh, job_single.sh)
scripts/probe_build.sh   o que a build do CP2K aceita
legacy/                  pipeline antigo (não usar; mantido como histórico)
```

Exemplo (bulk):

```bash
conda activate cp2k_env
python ~/work_cp2k/01_grid_convergence.py --phase cutoff --structure bulk.vasp \
    --basis DZVP-MOLOPT-SR-GTH --project X --values 300 400 500 600 700 800 --output-dir 01_grid_convergence
sbatch --array=0-5 ~/work_cp2k/jobs/job_scan_array.sh $PWD/01_grid_convergence/cutoff
python ~/work_cp2k/parsers/parse_grid.py --scan-dir 01_grid_convergence/cutoff
```
