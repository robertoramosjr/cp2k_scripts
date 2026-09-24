# CP2K_PIPELINE: especificação unificada do pipeline CP2K

> **Unificado em 2026-09-24 (máquina `access`)** a partir de duas fontes:
> (1) a spec original, escrita no `coaraci` em 2026-09-21/22 para o TiO2
> (`legacy/docs/pipeline_cp2k_coaraci.md`, cópia literal); (2) a reconstrução feita no `access` em
> 2026-09-24, quando o coaraci estava inacessível (`legacy/docs/CP2K_PIPELINE_ACCESS_RECONSTRUCTION_2026-09-24.md`).
> As PARTES 1 (1.1–1.11), 4 e 5 vêm literalmente da spec do coaraci. 1.12–1.19, 2, 3 e 6 foram
> fundidas ou acrescentadas no access. Contexto operacional (máquinas, fatos medidos, HSE06) em
> [SECOND_BRAIN.md](SECOND_BRAIN.md).

Escopo original: TiO2 célula unitária (≤10 átomos), PBE → estrutura eletrônica. Estendido em
2026-09-24 para superfícies (slabs, energia de superfície, Wulff), com ZrO2 c/t/m.

---

## PARTE 1: achados (o que estava errado e por quê)

Os itens 1.1–1.11 foram verificados na documentação oficial do CP2K (manual.cp2k.org,
cp2k.org/exercises) e nos papers do pacote (arXiv:2508.15559, arXiv:2003.03868). Os itens
1.12–1.19 foram verificados na build `cp2k/2026.1` do access (schema `--xml` e runs reais).

### 1.1 PDOS + k-points é uma combinação INVÁLIDA (causa raiz do OOM)

> "At present, it is not possible to get the projected density of states when
> doing a K-Point calculation." — cp2k.org/exercises:2017_uzh_cmest:pdos

> "cp2k now for pdos only samples the Gamma point." — lista cp2k

Confirmado também pela issue #4854 do repositório, que lista "PDOS" entre as
funcionalidades **ainda solicitadas** para k-points.

**Implicação:** todo o passo de "densificar o k-mesh 2× para o DOS" (hábito
correto no VASP) pede ao CP2K algo que ele não implementa. Não é ineficiente —
é inválido.

**Rota correta:** PDOS em Gamma-only sobre supercélula, com OT.
A amostragem implícita da zona de Brillouin vem do tamanho da supercélula.

### 1.2 ADDED_MOS superdimensionado (mecanismo direto do estouro de memória)

O pipeline antigo calculava `ADDED_MOS` = nº de MOs ocupados (imitando
`NBANDS=NELECT` do VASP). Os exemplos oficiais do CP2K usam valores de
**uma ordem de grandeza menor**:

| Fonte | Sistema | ADDED_MOS |
|---|---|---|
| cp2k.org/exercises:common:bs | WO3 (bandas, k-points) | 2 |
| arXiv:2508.15559 (paper oficial) | grafeno (bandas) | 2 |
| pipeline antigo (errado) | TiO2 6 átomos | 18 |

Por que explode: para **cada** k-point é resolvido um problema de autovalor
complexo Hermitiano generalizado (K^k C^k = S^k C^k ε^k). Com malha 12×12×12
(1728 k-points) × funções de onda complexas × contagem de bandas dobrada, o
custo de memória multiplica em três eixos simultaneamente.

**Regra nova:** ADDED_MOS entre 2 e 10 para este tamanho de sistema.

### 1.3 Keywords obrigatórias/recomendadas que faltavam

| Keyword | Seção | Por quê |
|---|---|---|
| `EXTRAPOLATION USE_GUESS` | `&QS` | Marcado como **required for K-Point sampling** no exercício oficial |
| `&SMEAR` (FERMI_DIRAC, 300 K) | `&SCF` | Sem ele o CP2K **não imprime a energia de Fermi** em cálculo com k-points |
| `PARALLEL_GROUP_SIZE` | `&KPOINTS` | Controle de paralelização/memória de k-points |
| `EPS_ADAPT 0.01` | `&DIAGONALIZATION` | Usado no exemplo oficial de bandas |
| `CELL_REF` | `&SUBSYS%CELL` | Ver 1.4 |

### 1.4 CELL_OPT sem CELL_REF gera descontinuidades artificiais

> "A variable-cell optimization may witness significant variations of the cell,
> which will also affect the number of grid points and introduce artificial
> discontinuities to the energy, force, stress and other properties even if the
> change in structure is small. It is recommended for better stability to set up
> a reference cell in CELL/CELL_REF, with the dimension of the reference cell
> larger than the expected upper bound […] which will be kept constant so that
> the number of grid points is held fixed."
> — manual, *Geometry and cell optimization*

O pipeline antigo não tinha `CELL_REF`. Todo CELL_OPT rodado até agora sofreu
esse problema.

### 1.5 A relaxação NÃO precisa de dois passos

`RUN_TYPE CELL_OPT` no CP2K já é definido como:

> "CELL_OPT — Cell optimization. **Both cell vectors and atomic positions are
> optimised.**" — manual, GLOBAL%RUN_TYPE

A divisão stress (ISIF=3) → Etot (ISIF=2) é uma convenção do VASP, não uma
necessidade do CP2K. Um único `CELL_OPT` bem convergido substitui os dois passos.

**Mapeamento correto:**

| VASP | CP2K | Relaxa |
|---|---|---|
| ISIF=3 | `CELL_OPT` | célula + íons (**suficiente sozinho**) |
| ISIF=2 | `GEO_OPT` | só íons (usar só se a célula for fixa por design) |

### 1.6 A convergência de CUTOFF estava sendo feita pelo critério errado

O método do próprio CP2K (manual, *How to Converge the CUTOFF and REL_CUTOFF*)
**não é** só olhar energia total vs. cutoff. É um procedimento 2D que também
analisa a **distribuição das Gaussianas pelos níveis do multigrid**:

1. `RUN_TYPE ENERGY` + `MAX_SCF 1` (sem SCF — muito mais barato que o que
   fazíamos com `ENERGY_FORCE` + SCF completo)
2. `PRINT_LEVEL MEDIUM` para imprimir quantas Gaussianas caem em cada grid
3. Varre CUTOFF com REL_CUTOFF fixo → escolhe pelo critério:

> "it is the lowest cutoff energy where the finest grid level is used, but at
> the same time with the majority of the Gaussians on the coarser grids"

4. **Depois** varre REL_CUTOFF com o CUTOFF escolhido fixo

Aviso importante do manual sobre a interação entre os dois:

> "simply increasing CUTOFF without increasing REL_CUTOFF may eventually lead to
> a slow convergence in energy, as more and more Gaussians get pushed to coarser
> grid levels, negating the increase in CUTOFF"

Ou seja: os dois parâmetros são acoplados. O pipeline antigo tratava REL_CUTOFF
como constante arbitrária (60) e nunca o convergiu.

### 1.7 O fator 1.5× de REL_CUTOFF para stress era invenção minha

Não tem respaldo na documentação. O que o manual recomenda para estabilidade de
CELL_OPT é `CELL_REF` (1.4), não margem no REL_CUTOFF. **Removido.**

### 1.8 HSE06 em bandas: limitação a verificar

> "K-point sampling only support GGA functionals. Band structure calculations
> using high-level electronic structure theory, such as hybrid functionals are
> not available in CP2K." — cp2k.org/exercises:common:bs (2023)

Documentação mais recente sugere que RI-HFX com k-points já existe ("HFX-RI with
k-Points […] includes a band-structure example"), mas **precisa ser testado na
build 2026.1 do cluster** antes de virar plano. O roadmap original que previa
bandas HSE06 direto pode não ser executável.

### 1.9 Interface nova de DOS: provavelmente indisponível

> "Most of the DOS/PDOS functionality described below, including the unified
> &DFT%PRINT%DOS interface and the &CURVE subsection […] is only available in
> **CP2K 2026.2 and later**."

Cluster tem **2026.1** → assumir indisponível; testar antes de depender.

### 1.10 BFGS vs LBFGS

> "BFGS — Most efficient minimizer, but only for 'small' systems, as it relies
> on diagonalization of a full Hessian matrix"

Para TiO2 unitário (≤10 átomos), BFGS é adequado. Para os sistemas de produção
(80-120 átomos, BaZrS3 etc.), trocar para LBFGS.

### 1.11 Bandas: a parte que estava certa

O uso de `&PRINT%BAND_STRUCTURE` com `&KPOINT_SET` está correto e é
computacionalmente barato:

> "any path of any density can be generated, since all relevant quantities are
> calculated in real-space" — arXiv:2508.15559

É o análogo funcional do `ICHARG=11` do VASP: uma SCF só, bandas obtidas por
transformada de Fourier das matrizes KS do espaço real. **Manter.**

### 1.12 (access) MPI: `srun` puro roda N cópias de 1 rank
O `cp2k/2026.1` do access usa o MPICH 4.3.2 do próprio toolchain. `srun` sem plugin PMI sobe N
cópias independentes (`Total number of message passing processes 1`), todas escrevendo no mesmo
`.out`; no TiO2_smoke do access foi isso que aconteceu, e os tempos de lá valem para 1 rank × 2
threads. **Lançar com o `mpiexec` do MPICH** (Hydra, SLURM-aware) e checar sempre essa linha. (O
mesmo sintoma apareceu no coaraci, onde a correção foi `srun --mpi=pmix`.)

### 1.13 (access) `&KPOINTS` em Γ aborta o OT
Escrever `&KPOINTS`, mesmo `MONKHORST-PACK 1 1 1`, põe o CP2K no caminho de k-points ("OT not
possible with kpoint calculations"). Em Γ-only a seção é omitida.

### 1.14 SCF não convergido aborta o run
`IGNORE_CONVERGENCE_FAILURE` é False por padrão (2025.1 e 2026.1). O método 2D (1.6), com
`MAX_SCF 1`, precisa ligá-lo; nenhuma outra etapa deve ligar.

### 1.15 `EXTERNAL_PRESSURE` default = 100 bar
No `&MOTION%CELL_OPT` (schema 2025.1 e 2026.1). Sem `EXTERNAL_PRESSURE [bar] 0` explícito, todo
CELL_OPT relaxa sob compressão.

### 1.16 O CP2K acrescenta ao `.out` existente
Parsers e jobs leem só o último segmento `PROGRAM STARTED AT`; ABORTs antigos já causaram
falso alarme. Os jobs renomeiam o `.out` antigo ao retomar.

### 1.17 O restart guarda `MAX_ITER` e `STEP_START_VAL`
Retomar depois de bater em MAX_ITER não faz nada até aumentar o MAX_ITER **dentro do
`.restart`**. O `.restart` também aninha `&CELL_REF` dentro de `&CELL`: um leitor ingênuo
("últimos A/B/C") devolve a célula 1.15×.

### 1.18 PDOS: formatos por versão
2026.1: só `&DFT%PRINT%PDOS` (o `&DFT%PRINT%DOS` só dá DOS total, `DELTA_E`). 2026.2:
`&DFT%PRINT%DOS` com `NLUMO` e `&PDOS` aninhado; o `&PRINT%PDOS` isolado deixa de existir e o
cabeçalho do `.pdos` muda (ver 5.4).

### 1.19 (access) Superfícies com pymatgen: três armadilhas
`Slab.is_polar()` é sempre False sem estados de oxidação; `SlabGenerator` e o construtor `Slab`
**giram a rede** por padrão (`reorient_lattice=True`), o que quebra a correspondência com os
vetores [uvw] do bulk; `get_tasker2_slabs()` falha em alguns planos Tasker III (t-ZrO2 (110)).

---

## PARTE 2: pipeline

```
[01] CUTOFF + REL_CUTOFF   01_grid_convergence.py   (método 2D, MAX_SCF 1)          ← "Etapa 0" do coaraci
        ↓
[02] k-mesh                02_kmesh_convergence.py  (CELL_OPT completo por malha)    ← "Etapa 1"
        ↓
[03] CELL_OPT de produção  03_cellopt.py            (passo único, célula + íons)     ← "Etapa 2"
        ↓
    ├── [04] bandas        04_bands.py              (k-mesh do 03 + BAND_STRUCTURE)  ← "Etapa 3a"
    ├── [05] PDOS          05_pdos.py               (supercélula, Γ, OT)             ← "Etapa 3b"
    └── [06] slabs         06_slab_cut.py → [07] 07_slab_opt.py → energia de superfície / Wulff  (PARTE 6)
Só no coaraci (não portados): 00_basis_probe.py, 06_bands_hybrid.py (entra aqui como 08_bands_hybrid.py).
```

`core/cp2k_blocks.py` é a fonte única dos blocos. Toda física tunável vem de `DEFAULTS` via
argparse, e `--basis` é obrigatório. Só os `jobs/*.sh` têm valores fixos. Os scripts rodam de
qualquer diretório e cada etapa grava `run_meta.json`/`scan_meta.json` para a próxima.

### Blocos comuns
- Γ-only: `&OT` (DIIS, FULL_SINGLE_INVERSE) + `&OUTER_SCF`, **sem** `&KPOINTS` (1.13).
- k ≠ Γ: `&DIAGONALIZATION` (STANDARD, `EPS_ADAPT 0.01`) + `&MIXING` (Broyden, α 0.4; 0.2 nos
  slabs) + `&SMEAR ON` (FERMI_DIRAC, 300 K) + `ADDED_MOS 10` (1.2, 1.3).
- `&QS`: `EPS_DEFAULT 1E-12`, `EXTRAPOLATION USE_GUESS`. `&KPOINTS PARALLEL_GROUP_SIZE -1`.
- `NGRIDS 5` (o coaraci usava 4; divergência mantida, registrada no SECOND_BRAIN).

### 01: CUTOFF e REL_CUTOFF (método oficial, 1.6)
`RUN_TYPE ENERGY`, `MAX_SCF 1` + `IGNORE_CONVERGENCE_FAILURE`, `PRINT_LEVEL MEDIUM`, Γ-only.
Fase `--phase cutoff`: CUTOFF 200/300→1000 Ry com REL_CUTOFF 60. Fase `--phase relcutoff`:
REL_CUTOFF 10→100 Ry com o CUTOFF escolhido. `parse_grid.py`: ΔE < limiar em meV/átomo, para o
ponto e para todos os maiores, **e** grid 1 usado com a maioria das Gaussianas nos grids grossos.
(No TiO2 do coaraci, o critério de 1e-8 Ha do manual era inatingível: o ruído é de ~4e-6 Ha.)

### 02: k-mesh
CELL_OPT completo por malha (`--k-list` ou `--densities`), critério duplo ΔE (meV/átomo) **e**
ΔV (%) contra a malha mais densa. Arrays usam `&SCF%PRINT%RESTART OFF` (o writer `.kp` crasha sob
I/O concorrente). `parse_kmesh.py` também imprime a densidade n_i·|a_i|, reusada por 04 e 07.

### 03: CELL_OPT de produção (passo único, 1.5)
`CELL_REF` = 1.15× os vetores (1.4), `STRESS_TENSOR ANALYTICAL`, `EXTERNAL_PRESSURE 0` (1.15),
BFGS até 10 átomos e LBFGS acima (1.10), `MAX_FORCE [eV*angstrom^-1] 0.010`,
`--keep-space-group` opcional, `&SCF%PRINT%RESTART ON`. `parse_cellopt.py` gera o CIF relaxado e o
`bulk.json` (E/f.u., k-density), que é a referência das energias de superfície.

### 04: bandas (1.11)
`RUN_TYPE ENERGY` na estrutura relaxada, SCF no k-mesh do 03 (convertido a densidade constante
se `kpath.prim` ≠ célula do 03), `ADDED_MOS 10` no SCF e no `&BAND_STRUCTURE`, caminho
`HighSymmKpath` com **`kpath.prim` no `&CELL`/`&COORD`**. `NPOINTS n` gera n+1 pontos por
segmento, um `KPOINT_SET` por ramo. `parse_bands.py` foi validado no c-ZrO2 real (2026.1): gap
indireto X→Γ de 3,30 eV.

### 05: PDOS (1.1)
Supercélula `--supercell NX NY NZ`, Γ-only, OT, `COMPONENTS`, `NLUMO -1`.
`--dos-interface legacy` (2025.1/2026.1) | `unified` (≥ 2026.2). O broadening fica em Python
(`parse_pdos.py --npoints 6000 --sigma`), e o parser lê os cabeçalhos dos dois formatos.

### Jobs
`jobs/common.sh` + `job_scan_array.sh` (01, 02, lotes de slabs; um subdiretório por tarefa,
`sort -V`) + `job_single.sh` (03, 04, 05, slab isolado). Padrão: `-N 1 --ntasks-per-node=14
--cpus-per-task=2 --mem-per-cpu=2194M`, `module cp2k/2026.1`, `CP2K_DATA_DIR="$CP2K_DATA"`,
`mpiexec -n $SLURM_NTASKS -bind-to core:2` (1.12), aviso se ranks ≠ `SLURM_NTASKS`. Restart-safe
via `<project>-1.restart`. Sucesso = `OPTIMIZATION COMPLETED` (GEO/CELL_OPT) ou `PROGRAM ENDED`
com energia (ENERGY); terminar não basta.

---

## PARTE 3: pendências (status em 2026-09-24)

1. `&DFT%PRINT%DOS` na 2026.1? **Respondido:** só DOS total; PDOS pela rota 3b (`legacy`).
   `scripts/probe_build.sh` confere isso em qualquer build.
2. RI-HFX com k-points na 2026.1? Sintaxe OK (`--check`). Falta o run real (`probe_build.sh --run`,
   só via sbatch). Na 2026.2 do coaraci funciona, mas o limite de memória torna inviável (SECOND_BRAIN §6).
3. Marcador de convergência: **respondido**, `OPTIMIZATION COMPLETED` (`GEOMETRY OPTIMIZATION COMPLETED`).
4. Ambiente conda: **respondido**, `cp2k_env` no access (pymatgen 2026.7.31).
5. Tamanho da supercélula do PDOS: em aberto (precisa de convergência própria).
6. Slabs: convergência do vácuo (**em fila**, array 5209911) e da espessura (em aberto);
   escolha de 2 planos cúbicos extras (com o usuário).
7. Detectar `[max_iter]` nos jobs do access (1.17); o coaraci já fazia isso.
8. Portar do coaraci: `00_basis_probe`, `08_bands_hybrid`, registro `BASES`, `run_smoke.py`,
   parsers `--json`, `&HF/&MEMORY` (5.3).

---

## PARTE 4: referências

- manual.cp2k.org/trunk/methods/dft/cutoff.html — convergência CUTOFF/REL_CUTOFF
- manual.cp2k.org/trunk/methods/optimization/geometry.html — CELL_REF, geo/cell opt
- manual.cp2k.org/trunk/methods/dft/k-points.html — k-points, wavefunctions complexas
- manual.cp2k.org/trunk/methods/electronic_structure/dos.html — DOS/PDOS, limitação de versão
- cp2k.org/exercises:common:bs — exemplo oficial de bandas (WO3)
- cp2k.org/exercises:2017_uzh_cmest:pdos — limitação PDOS + k-points
- github.com/cp2k/cp2k/issues/4854 — status de suporte a k-points
- arXiv:2508.15559 — *The CP2K Program Package Made Simple*
- arXiv:2003.03868 — *CP2K/Quickstep*
- pymatgen `SlabGenerator` / Tasker, *J. Phys. C* 12, 4977 (1979): classificação de superfícies iônicas (I/II/III).
- ZrO2 (planos de literatura): Christensen & Carter, PRB 58, 8050 (1998); Piskorz et al., JPCC 115, 24274 (2011) e 116, 19307 (2012).

---

## PARTE 5: adendo 2026-09-22 (coaraci): `&HF/&MEMORY`, `T_C_G_DATA` e a interface nova de PDOS

Motivado por um exemplo de HFX (PBE0, `&HF/&MEMORY MAX_MEMORY`, `T_C_G_DATA`,
`&INTERACTION_POTENTIAL TRUNCATED`) encontrado numa busca rápida, que levantou
a suspeita de estarmos "errando" algo no processo. Investigado e testado
contra o schema e o binário reais do CP2K 2026.2 instalado localmente.

### 5.1 O snippet encontrado tem um erro de sintaxe: `EPS_STORAGE_FILE` não existe

Reproduzido: um input com `&HF/&MEMORY/EPS_STORAGE_FILE` literal, como no
exemplo, faz o CP2K **abortar** ao ler o input:

```
[ABORT] found an unknown keyword EPS_STORAGE_FILE in section MEMORY
```

A keyword correta é `EPS_STORAGE` (alias `EPS_STORAGE_SCALING`), default 1.0
(escala de `EPS_SCHWARZ`; não precisa ser escrita explicitamente na maioria
dos casos). Corrigido e reverificado (`cp2k.ssmp --check` aceita).

### 5.2 `T_C_G_DATA`: não faltava — já resolve pelo diretório de dados

O default do CP2K para essa keyword já é o nome de arquivo `t_c_g.dat`,
resolvido pelo mesmo mecanismo de busca de `BASIS_SET_FILE_NAME` /
`POTENTIAL_FILE_NAME` (via `CP2K_DATA_DIR`). Confirmado: uma rodada real
RI-HFXk com `POTENTIAL_TYPE TRUNCATED`, sem cópia local do arquivo e sem
declarar `T_C_G_DATA` no input, terminou normalmente. **Não é preciso
declarar essa keyword.**

### 5.3 `&HF/&MEMORY` era um gap real — mas só no caminho Γ (HFX convencional), não no RI-HFXk de k-points

O exemplo do snippet usa HFX **convencional** (analítico, 4 centros,
`hfx_types.F`), o único caminho onde `&HF/&MEMORY` existe e faz sentido:
`MAX_MEMORY` define o cache de integrais por rank MPI (default do próprio
CP2K: **512 MiB**, baixo). O pipeline não escrevia essa seção no caminho
Γ-only híbrido (`build_xc_block`, usado por `05_pdos.py` com PBE0/HSE06 em
Γ) — corrigido: agora escreve `&HF/&MEMORY MAX_MEMORY <valor>` sempre nesse
caminho, default 4096 MiB/rank, ajustável via `--hfx-max-memory`.

**Isso NÃO é a causa do OOM de bandas híbridas com k-points** (Seção 1.8 e o
adendo de bandas híbridas, 2026-09-21): RI-HFXk usa `&HF/&RI` — uma seção
IRMÃ de `&MEMORY`, não filha — com seu próprio modelo de memória
(`KP_NGROUPS`, `KP_STACK_SIZE`, `MEMORY_CUT`, já implementados), porque RI
decompõe as integrais de 4 centros via base auxiliar e nunca monta/armazena
a lista de ERIs de 4 centros que `&MEMORY` controla. `&HF/&MEMORY` é
ignorado nesse caminho e não foi adicionado lá.

### 5.4 Achado à parte (verificando o item 1.9): `&PRINT%PDOS` mudou de lugar no CP2K 2026.2

Ao testar `05_pdos.py` de ponta a ponta no CP2K 2026.2 local, o `&DFT%PRINT%PDOS`
isolado (o único formato que o pipeline gerava) **deixou de existir** nessa
versão — `cp2k.ssmp --check` aborta com "unknown subsection PDOS of section
PRINT". A seção foi movida para dentro da interface nova de DOS unificada
(`&DFT%PRINT%DOS%PDOS`), com `NLUMO` subindo para o nível de `&DOS`. Esse é
exatamente o cenário que o item 1.9 pediu para testar antes de depender.

Corrigido com um parâmetro `--dos-interface {legacy,unified}` em
`05_pdos.py` (`legacy` continua default, é o formato verificado contra
rodadas reais em CP2K 2025.1; `unified` gera o formato novo). Também
descoberto no mesmo teste: o cabeçalho do arquivo `.pdos` em si mudou (nome
do kind e `E(Fermi)` passam a ficar em linhas separadas, e todas as
componentes até f aparecem mesmo sem elétrons f) — `parsers/parse_pdos.py`
foi ajustado para ler os dois formatos.

### 5.5 O que ainda não está resolvido

Nada disso muda o limite de memória do RI-HFXk com k-points (bandas
híbridas): esse continua bloqueado pelo tamanho da base RI estendida, não
por `&HF/&MEMORY` nem por `T_C_G_DATA`. Ver a memória de sessão
`cp2k_hybrid_bands_addendum.md` para os números medidos.

---

## PARTE 6: superfícies (access, 2026-09-24)

### 06_slab_cut.py
A partir do bulk convergido (o CONTCAR do VASP define o setting dos índices; em P2₁/c, (1̄11) ≠
(111)), para cada plano: `SlabGenerator` (`reorient_lattice=False`) com as terminações simples e
simetrizadas, mais as reconstruções Tasker III→II (a do pymatgen e `tasker3_to_2`: metade da camada
externa movida um período de empilhamento para o lado oposto). Ficam só os slabs
**estequiométricos, apolares (cargas formais, `--oxidation`) e simétricos**, sem duplicatas,
ordenados por ligações cátion–ânion quebradas por Å². Gera `slab_<hkl>_<n>/` com c ⟂ superfície
∥ z e `slab_meta.json` (hkl, n de f.u., área, [uvw] racionais do plano, autoteste do hkl pela
normal).

### 07_slab_opt.py
Reconstrói os vetores do plano com os mesmos [uvw] sobre a rede do 03 e escala a altura por
d_hkl(CP2K)/d_hkl(ref), o que remove o descasamento VASP→CP2K como deformação artificial.
`--vacuum` escolhe o vácuo; `--run-type GEO_OPT` (produção, célula fixa) ou `ENERGY` (testes).
`SURFACE_DIPOLE_CORRECTION`/`SURF_DIP_DIR Z` fica ligado por segurança (dipolo nulo nos slabs
simétricos); k-mesh `ceil(k_density/|a|) × ceil(k_density/|b|) × 1`; α de mistura 0.2.

### Energia de superfície e Wulff
`parse_surface_energy.py`: γ = (E_slab − n·E_bulk/f.u.)/(2A) em J/m², checando se cutoff, basis e
k-density batem com o `bulk.json`; a melhor terminação por plano entra no `WulffShape` (frações
de área = morfologia). `parse_vacuum.py`: Δγ(V) = [E(V) − E(V_max)]/(2A) por slab; o
vácuo recomendado é o menor com |Δγ| < 5 mJ/m² para ele e todos os maiores.

### Vácuo (por que é necessário e quanto custa)
A célula é periódica em XYZ: sem vácuo o slab vira bulk com defeito de empilhamento, em qualquer
código. No GPW o vácuo não acrescenta funções de base (os orbitais são Gaussianas centradas nos
átomos) e só aumenta a grade da densidade (FFT/Hartree/XC). A parte cara (diagonalização em
k-points) não depende dele. `PERIODIC XY` com solver 2D não resolve: o MT exige célula ≥ 2× a
extensão da densidade, e o wavelet exige célula ortorrômbica.
