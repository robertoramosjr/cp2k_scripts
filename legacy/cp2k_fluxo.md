# Fluxo CP2K -- TiO2 Anatase (convergência de parâmetros)

Manual de uso dos scripts em `~/work_cp2k/scripts/`. Pensado tanto para quem já rodou
esse pipeline quanto para outra pessoa pegando o projeto pela primeira vez.

**Ordem de dependência:** `CUTOFF -> k-mesh (stress, rigoroso) -> extrair estrutura relaxada -> Etot (produção) -> DOS + Bandas`

Cada passo alimenta o próximo: o `CUTOFF` recomendado do passo 1 é fixado no passo 2;
a malha k recomendada e a estrutura relaxada do passo 2 alimentam o passo 4; a
estrutura e a malha k do Etot alimentam o passo 5.

---

## 0. Antes de começar (setup do ambiente)

```bash
module purge
module load miniconda
module load intel/mpi/2017
module load cp2k/2026.1
```

O módulo `cp2k/2026.1` exporta `CP2K_DATA` (nome não-padrão), mas o binário
(`cp2k.psmp`) só lê `CP2K_DATA_DIR`. **Sempre exportar isso no job**, ou o CP2K
cai no default `./data` (relativo ao diretório de execução) e aborta com
`The specified OLD file <BASIS_MOLOPT> cannot be opened`:

```bash
export CP2K_DATA_DIR="$CP2K_DATA"
```

Já embutido em todos os scripts de job listados na seção 7.

**Materials Project:** defina `MP_API_KEY` no ambiente (ou passe `--api-key`) antes de
usar `--mp-id` em qualquer gerador. Se `python3 -c "from mp_api.client import MPRester"`
der traceback, é conflito de versão `mp-api`/`pymatgen`/`emmet-core` -- considere um
`cp2k_env` conda separado do `vasp_env`.

---

## 1. Convergência de CUTOFF (barata, ponto único)

Equivalente a um scan de `ENCUT` no VASP: varia o `MGRID CUTOFF`, k-mesh fixo
(Gamma por padrão) e barato -- `RUN_TYPE ENERGY_FORCE`, sem otimizar nada.

```bash
python ~/work_cp2k/scripts/cutoff_convergence.py --mp-id mp-390 --project-name TiO2_anatase \
    --cutoff-min 200 --cutoff-max 1000 --cutoff-step 100
```

**Rodar -- opção A, loop linear** (1 `sbatch`, todos os pontos em série):
```bash
sbatch ~/work_cp2k/scripts/job_cutoff_convergence.sh $(pwd)/cutoff_convergence
```

**Rodar -- opção B, job array** (pontos em paralelo, melhor throughput):
```bash
N=$(find cutoff_convergence -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=0-$((N-1)) ~/work_cp2k/scripts/job_cutoff_convergence_array.sh $(pwd)/cutoff_convergence

# resubmissão de 1 índice específico que falhou/não terminou:
sbatch --array=5 ~/work_cp2k/scripts/job_cutoff_convergence_array.sh $(pwd)/cutoff_convergence
```
Status por ponto fica em `cutoff_convergence/<dir>/status.txt` (não `calculation.dat` --
evita race condition entre tarefas concorrentes escrevendo no mesmo arquivo).

**Parse:**
```bash
python ~/work_cp2k/scripts/parse_cutoff_convergence.py --scan-dir cutoff_convergence --n-atoms <N_ATOMS>
```
-> anota o **CUTOFF recomendado** (primeiro valor dentro do threshold de meV/átomo
em relação ao cutoff mais alto testado).

> `--n-atoms` tem que ser o número real de átomos da célula que você está rodando
> (o próprio gerador imprime isso: `Atoms in cell (needed later for meV/atom
> normalization): N`) -- usar um valor errado (ex: de uma célula primitiva menor)
> distorce o `dE (meV/átomo)` por um fator igual à razão entre os dois.

---

## 2. Convergência de k-mesh via STRESS completo (metodologia rigorosa)

Replica o job de "stress" do VASP: em vez de um ponto único, cada malha k testada
roda um `CELL_OPT` completo (célula + íons relaxando juntos, `ISIF=3` equivalente),
usando o `CUTOFF` já convergido no passo 1. Mais caro, mas testa a grandeza física
que importa (energia relaxada + volume relaxado) em vez de energia numa geometria
ainda não relaxada.

```bash
python ~/work_cp2k/scripts/kmesh_stress_convergence.py --mp-id mp-390 --project-name TiO2_anatase \
    --cutoff 600 --k-min 2 --k-max 10 --k-step 2

# ou uma lista explícita de malhas não-isotrópicas:
python ~/work_cp2k/scripts/kmesh_stress_convergence.py --mp-id mp-390 --project-name TiO2_anatase \
    --cutoff 600 --k-list "2,1,1;6,2,2;6,4,4;9,4,4;12,5,5"
```

**Rodar -- opção A, loop linear:**
```bash
sbatch ~/work_cp2k/scripts/job_kmesh_stress_convergence.sh $(pwd)/kmesh_stress_convergence
# Se estourar walltime: resubmete o MESMO sbatch. O job detecta o .restart de cada
# diretório automaticamente e continua de onde parou -- diferente do VASP, aqui não
# precisa copiar CONTCAR->POSCAR manualmente pra RETOMAR o mesmo cálculo: o .restart
# do CP2K já é um input completo por si só (isso é só pra resume; ver seção 3 para
# a extração que alimenta o PRÓXIMO passo do pipeline, que aí sim é como cp CONTCAR POSCAR).
```

**Rodar -- opção B, job array:**
```bash
N=$(find kmesh_stress_convergence -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=0-$((N-1)) ~/work_cp2k/scripts/job_kmesh_stress_convergence_array.sh \
    $(pwd)/kmesh_stress_convergence

# resubmissão de 1 índice específico (mesma lógica de restart por diretório):
sbatch --array=3 ~/work_cp2k/scripts/job_kmesh_stress_convergence_array.sh \
    $(pwd)/kmesh_stress_convergence
```
Pesa ainda mais aqui a escolha loop vs. array do que no passo 1, já que cada ponto
é um `CELL_OPT` completo (bem mais caro que um ponto único) -- ver seção 7.

**Convergência real vs. "terminou sem crashar":** os dois jobs acima checam
`RUN_TYPE` do próprio `.inp` antes de decidir o critério de sucesso -- banner de
otimização (`OPTIMIZATION COMPLETED`) para `CELL_OPT`/`GEO_OPT`, ou
`outer SCF loop converged` para `ENERGY_FORCE`/`ENERGY`. Isso importa porque o CP2K
**não aborta** quando o SCF não converge dentro do `MAX_SCF`/`OUTER_SCF` -- ele só
segue em frente e imprime a energia mesmo assim, sem erro nenhum visível no
`.inp.out`. Um `PROGRAM ENDED AT` sozinho NÃO significa que o resultado é confiável.

**Parse:**
```bash
python ~/work_cp2k/scripts/parse_kmesh_stress_convergence.py --scan-dir kmesh_stress_convergence \
    --n-atoms <N_ATOMS> --energy-threshold-mev-atom 1.0 --volume-threshold-pct 0.5
```
-> anota a **malha k recomendada** (ex: `6x6x6`), escolhida cruzando ΔE (meV/átomo)
**e** ΔVolume (%) contra a malha mais densa testada.

Detalhes técnicos aplicados automaticamente nesse passo (via `input_generator.py`):
- `STRESS_TENSOR ANALYTICAL` no `&FORCE_EVAL` -- necessário pro CP2K calcular o
  tensor de stress e variar a célula; sem isso o `CELL_OPT` aborta.
- `REL_CUTOFF` ampliado por `--stress-rel-cutoff-factor` (default 1.5x) -- análogo
  do CP2K/GPW ao hábito de dar margem extra na etapa de volume no VASP, mas aplicado
  ao `REL_CUTOFF` (completude do multigrid), não ao `CUTOFF` puro, já que a base
  Gaussiana não depende do volume da célula como ondas planas dependem.
- Malhas k não-Gamma usam `&DIAGONALIZATION` + `&MIXING`, nunca `&OT` -- o
  minimizador OT do CP2K não suporta k-points (aborta com `OT not possible with
  kpoint calculations`), mesmo numa malha `1x1x1` explícita. Detectado automaticamente
  pelos geradores (`is_gamma_only`); só `1x1x1` de fato usa `&OT`.
- `&MOTION &PRINT &RESTART` com `EACH CELL_OPT 1` -- garante checkpoint a cada passo.
- Restart de wavefunction (`&SCF/&PRINT/&RESTART`) desativado automaticamente em
  qualquer malha não-Gamma, mesmo em `CELL_OPT`/`GEO_OPT` -- ver pegadinha na seção 8
  sobre o crash `SIGABRT` em `write_kpoints_restart`.

---

## 3. Extrair a estrutura relaxada do k-mesh escolhido

Isto é o `cp CONTCAR POSCAR` do fluxo: pega o resultado final de um cálculo
(`CELL_OPT` da malha vencedora) e usa como ponto de partida de um cálculo
**diferente** (o `GEO_OPT` de produção do passo 4). Diferente do resume dentro do
mesmo cálculo (passo 2), aqui é necessário porque o `.restart` do CP2K é um
checkpoint mais "gordo" que o `CONTCAR` (carrega estado de SCF, contadores de
passo etc.) -- o script extrai só a parte equivalente ao `CONTCAR` (célula +
coordenadas finais) e reempacota como CIF limpo.

```bash
python ~/work_cp2k/scripts/extract_relaxed_structure.py \
    --restart kmesh_stress_convergence/k6-6-6/TiO2_anatase_k6-6-6-1.restart \
    --output TiO2_anatase_relaxed.cif
```

---

## 4. Etot -- relaxação iônica de produção (célula fixa, `ISIF=2` equivalente)

```bash
python ~/work_cp2k/scripts/input_generator.py --cif TiO2_anatase_relaxed.cif \
    --project-name TiO2_etot --run-type GEO_OPT --cutoff 600 --kpoints 6 6 6
```
`RUN_TYPE GEO_OPT` nunca mexe na célula -- só existe `&MOTION/&GEO_OPT`, sem
`&CELL_OPT` nem `STRESS_TENSOR`. Célula fixa no volume/forma vindos do passo 3,
só as posições atômicas relaxam. Sem margem extra de `REL_CUTOFF` aqui -- célula já
fixa, não precisa da precisão extra que o tensor de stress exige.

**Rodar:**
```bash
sbatch ~/work_cp2k/scripts/job_etot.sh $(pwd)/etot_dir
```
Job único (sem array -- é uma relaxação só), com protocolo de restart igual aos
outros (`cp CONTCAR POSCAR` automático via `.restart` se estourar walltime).

---

## 5. DOS + Bandas (PBE, base para reaproveitar via restart em HSE06 depois)

Uma única run de CP2K calcula PDOS e bandas juntas (mesma SCF autoconsistente;
bandas são interpoladas por Fourier das matrizes de Kohn-Sham no espaço real,
não precisam de segunda SCF -- equivalente ao `ICHARG=11` do VASP).

### 5.1 Puxar os arquivos do Etot pra pasta nova (mantém histórico separado)

```bash
mkdir dos_bands_TiO2 && cd dos_bands_TiO2

cp /caminho/etot_run/TiO2_etot-1.restart .
cp /caminho/etot_run/TiO2_etot-RESTART.wfn .

python ~/scripts/extract_relaxed_structure.py --restart TiO2_etot-1.restart \
    --output TiO2_etot_relaxed.cif
```

### 5.2 Rodar o pipeline completo (1 job, 1 alocação de fila)

Edita o bloco `USER CONFIG` no topo do `job_dos_bands.sh` antes de submeter:

```bash
RELAXED_CIF="TiO2_etot_relaxed.cif"    # confira maiúsculas/minúsculas contra o
                                        # arquivo real -- Linux é case-sensitive
                                        # (já pegou esse erro uma vez neste projeto)
PROJECT_NAME="TiO2"
CUTOFF=<MESMO CUTOFF DO ETOT>          # não muda entre Etot e DOS/bandas
KPOINTS="<MESMO K-MESH DO ETOT>"       # ex: "6 6 6" -- é o BASELINE, não o mesh final
VERIFY_FUNCTIONAL="PBE"                # Etapa 1 (verify-relax) -- deve bater com o
                                        # funcional que RELAXOU o RELAXED_CIF (o Etot),
                                        # não necessariamente igual ao FUNCTIONAL abaixo
FUNCTIONAL="HSE06"                     # Etapas 3-4 (a run real de DOS/bandas)
DOS_KPOINTS_MULTIPLIER=2.0             # k-mesh real da run = KPOINTS x este fator
WFN_RESTART="TiO2_etot-RESTART.wfn"    # "" para pular SCF_GUESS RESTART
```

```bash
sbatch ~/scripts/job_dos_bands.sh $(pwd)
```

O job faz sozinho, dentro da mesma alocação:
1. **verify-relax**: `GEO_OPT` curto (funcional `VERIFY_FUNCTIONAL`) na célula
   padronizada pelo pymatgen (confirma que a padronização/simetrização não desviou
   do mínimo relaxado real do Etot).
2. Extrai a estrutura confirmada (`extract_relaxed_structure.py`).
3. Gera o `.inp` de DOS/bandas (`dos_bands_generator.py`, funcional `FUNCTIONAL`), com:
   - `ADDED_MOS` automático = elétrons de valência / 2 (equivalente ao
     `NBANDS=NELECT` do VASP).
   - k-mesh real = `KPOINTS x DOS_KPOINTS_MULTIPLIER` (arredondado) --
     PDOS precisa de malha mais densa que a relaxação, igual no VASP.
4. Roda o DOS/bandas.

`VERIFY_FUNCTIONAL` e `FUNCTIONAL` são propositalmente separados: a etapa 1 só
precisa confirmar o mínimo na MESMA superfície de energia que gerou `RELAXED_CIF`
(normalmente PBE, o funcional usado no Etot) -- rodar essa confirmação em HSE06
seria mais caro e não faz sentido fisicamente, já que a estrutura não foi relaxada
nesse nível de teoria. `FUNCTIONAL` (etapas 3-4) é o nível de teoria de interesse
pra DOS/bandas em si, que pode ser diferente (ex: HSE06 pra corrigir o gap
subestimado do PBE).

Restart-safe: se estourar walltime na etapa 1, resubmete o MESMO `sbatch` --
ele detecta o progresso e retoma sozinho. A etapa 4 (SCF de ponto único) checa
tanto `PROGRAM ENDED AT` quanto `outer SCF loop converged` antes de marcar
`OK` -- ver pegadinha na seção 8 sobre "terminou" vs. "convergiu de verdade".

### 5.3 Reaproveitando pra HSE06 depois

O `.wfn` do PBE pode alimentar `--wfn-restart` mesmo trocando `FUNCTIONAL` pra
`HSE06` -- a documentação do CP2K recomenda isso explicitamente pra acelerar a
SCF de híbridos (a base primária Gaussiana precisa ser a mesma; a base auxiliar
do ADMM é adicional, não quebra o restart).

---

## 6. Tabela de equivalência VASP -> CP2K

| VASP (INCAR) | CP2K | Onde aparece |
|---|---|---|
| `ISIF=3` (célula + íons juntos) | `RUN_TYPE CELL_OPT` | Passo 2 |
| `ISIF=2` (só íons, célula fixa) | `RUN_TYPE GEO_OPT` (sem `&CELL_OPT`) | Passo 4 |
| `ICHARG=11` (DOS/bandas na densidade convergida) | `RUN_TYPE ENERGY` + `&PRINT/&PDOS` + `&PRINT/&BAND_STRUCTURE` | Passo 5 |
| `NBANDS=NELECT` | `ADDED_MOS` = elétrons de valência / 2 (auto) | Passo 5 |
| `CONTCAR` | `<project>-1.restart` | Escrito por `&MOTION &PRINT &RESTART` |
| `cp CONTCAR POSCAR` (próximo passo) | `extract_relaxed_structure.py` | Passo 3, e antes do passo 5 |
| `WAVECAR` (seed pra próxima SCF) | `<project>-RESTART.wfn` | `SCF_GUESS RESTART` via `--wfn-restart` |
| `EDIFF = 1.0E-6` | `EPS_SCF 1.0E-6` | `&SCF` -- mesma ordem de grandeza, definição matemática difere levemente |
| `NELM = 1500` | `MAX_SCF 300` (dentro de `&OT`/`&DIAGONALIZATION`) | `&SCF` -- **não** é cópia 1:1; OT/DIIS do CP2K converge bem mais rápido que o Davidson do VASP, 300 é margem generosa |
| `AMIX = 0.4` | `&MIXING ALPHA 0.4` | Só se aplica em malhas k não-Gamma (`&DIAGONALIZATION`); Gamma usa `&OT`, que não tem mixing |
| `NSW = 800` | `MAX_ITER 500` | `&GEO_OPT`/`&CELL_OPT` -- mesma lógica de margem generosa (BFGS converge bem antes disso) |
| `EDIFFG = -0.01` | `MAX_FORCE [eV*angstrom^-1] 0.010` | `&GEO_OPT`/`&CELL_OPT` -- unidade nativa do CP2K, sem conversão manual |
| ENCUT scan | `cutoff_convergence.py` | Passo 1 |
| KPOINTS scan (stress) | `kmesh_stress_convergence.py` | Passo 2 |

Todos os valores da coluna CP2K acima (exceto DOS/bandas) são **configuráveis via
argparse** em `input_generator.py`, `cutoff_convergence.py` e
`kmesh_stress_convergence.py` -- os valores da tabela são os defaults, não
hardcoded. Flags disponíveis: `--eps-scf`, `--max-scf`, `--outer-max-scf`,
`--mixing-alpha`, `--max-iter`, `--max-force` (os dois últimos não existem em
`cutoff_convergence.py`, que é `ENERGY_FORCE` e não roda otimizador).

---

## 7. Referência dos scripts

### Geradores de `.inp`

| Script | Papel |
|---|---|
| `input_generator.py` | Motor compartilhado (biblioteca, `import input_generator as ig`) usado pelos outros geradores; também standalone para o Passo 4 (produção) ou qualquer `.inp` avulso fora do pipeline de convergência. |
| `cutoff_convergence.py` | Passo 1 -- gera N `.inp` `ENERGY_FORCE` variando `CUTOFF`. |
| `kmesh_stress_convergence.py` | Passo 2 -- gera N `.inp` `CELL_OPT` variando a malha k. |
| `extract_relaxed_structure.py` | Passos 3 e 5.1 -- extrai célula+coordenadas finais de um `.restart` para CIF. |
| `dos_bands_generator.py` | Passo 5 -- gera o `.inp` de verify-relax e o `.inp` combinado de PDOS + bandas (caminho de k-pontos de alta simetria via `pymatgen.HighSymmKpath`). |

### Jobs SLURM

| Script | Modo | Uso |
|---|---|---|
| `job_cutoff_convergence.sh` | Loop linear | 1 `sbatch`, todos os pontos em série |
| `job_cutoff_convergence_array.sh` | Job array | 1 `sbatch --array=0-N`, pontos em paralelo |
| `job_kmesh_stress_convergence.sh` | Loop linear | idem, pra scan de k-mesh (`CELL_OPT` completo) |
| `job_kmesh_stress_convergence_array.sh` | Job array | idem, versão array |
| `job_etot.sh` | Job único | relaxação de produção (Passo 4), sem array -- é uma run só |
| `job_dos_bands.sh` | Job único, 2 etapas de CP2K | Passo 5 completo (verify-relax + DOS/bandas) numa alocação só |

`-N 1 --ntasks-per-node=8 --cpus-per-task=2` (ou `-N 1 -n 8 -c 2`) é o ponto de
partida validado pros arrays, dimensionado pro sistema pequeno atual (anatase,
6-24 átomos dependendo da célula). **Quando usar qual:** array compensa pra
sistemas pequenos onde um cálculo individual não aproveita um nó inteiro -- melhor
rodar vários pontos ao mesmo tempo. Pra sistemas de produção grandes (~80-120
átomos), faça um teste de scaling (8/16/28/56 cores) primeiro; se o cálculo
individual escalar bem até perto do nó inteiro, o array ainda vale (é só um
mecanismo de agendamento), mas o `--ntasks-per-node`/`--cpus-per-task` de cada
tarefa deve subir proporcionalmente -- não faz sentido fatiar um sistema grande em
pedaços pequenos.

### Parsers

| Script | Lê a saída de | Critério de recomendação |
|---|---|---|
| `parse_cutoff_convergence.py` | Passo 1 | ΔE (meV/átomo) vs. cutoff mais alto testado |
| `parse_kmesh_stress_convergence.py` | Passo 2 | ΔE (meV/átomo) **e** ΔVolume (%) vs. malha mais densa testada |
| `parse_dos_bands.py` | Passo 5 | **Ainda não construído** -- ver pendências |

Os dois já existentes exigem que o SCF (e, no caso do passo 2, também o
otimizador) tenha convergido de verdade -- pontos sem convergência real são
pulados com `WARNING`, não incluídos silenciosamente na curva.

### Descontinuados (não usar pra k-mesh)

`kpoint_convergence.py` / `parse_kpoint_convergence.py` faziam a versão **barata**
(ponto único, `ENERGY_FORCE`, sem otimizar célula) da convergência de k-mesh --
substituídos por `kmesh_stress_convergence.py` / `parse_kmesh_stress_convergence.py`
por decisão de replicar a metodologia rigorosa do VASP (`CELL_OPT` completo). Não
estão mais neste diretório. Dados legados dessa rota antiga (se existirem) ficam
em diretórios como `kpoint_convergence/` -- pontos `ENERGY_FORCE`, sem relaxação
de célula, úteis só como pré-checagem rápida antes do scan completo de stress.

---

## 8. Pegadinhas já encontradas nesse projeto

- **`CP2K_DATA_DIR` não é setada pelo módulo** -- ver seção 0. Sintoma: abort
  silencioso em `read_qs_kind` procurando `BASIS_MOLOPT` em `./data`.
- **`&OT` nunca funciona com `&KPOINTS` explícito**, mesmo `1x1x1` -- os geradores
  já lidam com isso automaticamente (só omitem `&KPOINTS`/usam `&OT` quando o mesh
  é literalmente `1,1,1`; qualquer outro valor usa `&DIAGONALIZATION`+`&MIXING`).
- **`PROGRAM ENDED AT` != convergência real.** CP2K não aborta quando o SCF ou o
  otimizador não convergem -- os jobs em `job_*_array.sh`/`job_*.sh` e os parsers
  checam explicitamente `outer SCF loop converged` (ponto único) ou o banner de
  otimização (`CELL_OPT`/`GEO_OPT`), não só o fim do programa.
- **`--n-atoms` errado distorce meV/átomo silenciosamente** -- sempre conferir
  contra o que o próprio gerador imprime (`Atoms in cell...`), não assumir.
- **`sbatch --mail-type` precisa de `BEGIN` explícito** -- sem isso, nenhum e-mail
  de início é mandado (só fim/falha). Já corrigido em todos os jobs desta lista,
  mas vale checar em scripts novos copiados a partir de versões antigas.
- **Diretório de scan errado (`SCAN_DIR`) faz o loop iterar sobre a pasta errada**
  sem erro óbvio -- se `SCAN_DIR` apontar um nível acima do esperado, o loop trata
  a pasta inteira como um único "ponto", não acha `.inp` dentro dela, e sai sem
  rodar nada. Sempre apontar `SCAN_DIR` direto pra pasta que contém as subpastas
  `cutoffNNN/`/`kN-N-N/`, não a pasta pai.
- **Crash `SIGABRT` em `write_kpoints_restart` (k-mesh não-Gamma + `CELL_OPT`/
  `GEO_OPT`)** -- o restart de wavefunction (mesma keyword `&SCF/&PRINT/&RESTART`
  usada pro caso Gamma) escreve um `<project>-RESTART.kp` binário grande; sob
  concorrência de I/O (vários array tasks escrevendo ao mesmo tempo no mesmo
  filesystem) o CP2K crasha com `SIGABRT` dentro dessa rotina -- sem `[ABORT]`,
  sem nada visível no `.inp.out` (só aparece no `.error` do SLURM). Mitigado
  desativando esse restart também para qualquer malha não-Gamma, não só pra runs
  de ponto único. O checkpoint de geometria (`&MOTION/&PRINT/&RESTART`) não é
  afetado -- resume de walltime continua funcionando, só perde o "chute inicial"
  de densidade eletrônica ao retomar (recomeça de `SCF_GUESS ATOMIC`).
- **`sort` puro ordena diretórios errado em jobs array** -- `k12-3-3` vem antes de
  `k2-1-1` porque `sort` compara caractere por caractere (`'1' < '2'`), não como
  número. Os scripts de array usam `sort -V` (natural/version sort) pra dar a
  ordem numérica esperada.
- **Nunca edite um script de array enquanto instâncias dele ainda estão
  rodando/pendentes na mesma submissão** -- cada tarefa pendente lê o script do
  disco de novo quando começa (não herda o estado do momento do `sbatch`). Editar
  o `sort`/mapeamento de diretórios no meio de um array já em execução pode fazer
  duas tarefas caírem no mesmo diretório ao mesmo tempo. Sempre `squeue -u
  $USER` antes de editar um script associado a um job ativo; se precisar editar,
  cancele as tarefas ainda pendentes primeiro (as que já rodam não são afetadas,
  já carregaram o script em memória).
- **Essa build do CP2K 2026.1 não tem a subseção genérica `&XC_FUNCTIONAL/
  &LIBXC`** -- funcionais híbridos do libxc (ex: HSE06) precisam ser ativados
  como sua própria subseção direta (`&HYB_GGA_XC_HSE06 &END HYB_GGA_XC_HSE06`),
  não com `&LIBXC / FUNCTIONAL XC_HYB_GGA_XC_HSE06`. Sintoma:
  `[ABORT] unknown subsection LIBXC of section XC_FUNCTIONAL`. Pra descobrir a
  estrutura real de qualquer seção nessa build, gere o schema completo com
  `cp2k.psmp --xml` (escreve `cp2k_input.xml` no diretório atual) em vez de
  confiar em exemplos de versões antigas do CP2K.
- **Nomes de arquivo são case-sensitive no Linux** -- `TiO2_relaxed.cif` e
  `tio2_relaxed.cif` são arquivos diferentes. Sempre conferir a config
  (`RELAXED_CIF=...` etc.) contra `ls` do que realmente existe na pasta, não
  contra o que "deveria" estar lá.
- **Um `--functional` esquecido cai no default errado, sem avisar** -- a etapa
  1 (`verify-relax`) do `job_dos_bands.sh` chamava `dos_bands_generator.py` sem
  passar `--functional`, então rodava sempre com o default do gerador (`PBE`),
  mesmo com `FUNCTIONAL="HSE06"` configurado no job -- coincidentemente
  inofensivo aqui (verify-relax *devia* usar o mesmo nível de teoria do Etot,
  não o da run de DOS/bandas), mas implícito e frágil: se o default do gerador
  mudasse, o comportamento do job mudaria junto, sem nenhuma linha editada nele.
  Corrigido tornando isso explícito com duas variáveis separadas
  (`VERIFY_FUNCTIONAL` / `FUNCTIONAL`, seção 5.2) -- desconfie de qualquer
  script que dependa do default silencioso de outro script pra se comportar
  certo.

- **(2026-09-24) `srun` lançava N cópias seriais do mesmo cálculo.** O `cp2k/2026.1` é
  linkado contra o MPICH 4.3.2 do próprio toolchain, não contra Intel MPI; com `srun` sem
  plugin `--mpi` cada rank subia como um job de 1 processo (`Total number of message passing
  processes 1`), todos escrevendo no mesmo `.out`. Todos os runs antigos (TiO2_smoke)
  rodaram assim: tempos = 1 rank x 2 threads. Corrigido em todos os `job_*.sh`:
  `mpiexec -n $SLURM_NTASKS` do MPICH, sem `module load intel/mpi/2017`, e
  `check_mpi_ranks` avisa no `.error` se o nº de ranks do CP2K não bater com o SLURM.
- **(2026-09-24) `&KPOINTS` com malha 1x1x1 + `&OT` aborta** ("OT not possible with kpoint
  calculations"): `render_input` escrevia a seção sempre. Agora `build_kpoints_block()` omite
  a seção em Γ (vale também para `cutoff_convergence.py`).
- **(2026-09-24) `cutoff_convergence.py`/`kmesh_stress_convergence.py` passavam
  `mixing_alpha` para `build_scf_block`, que não aceitava o argumento** (drift entre
  scripts): TypeError na geração. Agora `--mixing-alpha` existe em `input_generator.py`.

---

## 9. Pendências em aberto

- `parse_dos_bands.py` (broadening Gaussiano tipo NEDOS=6000 sobre os `.pdos`
  crus, plot de bandas a partir do `.bs`) -- ainda não construído.
- Confirmar se a interface nativa `&DFT%PRINT%DOS%CURVE` (broadening no
  próprio CP2K, via `DELTA_E`) está disponível na build 2026.1 real do
  cluster antes de tentar usá-la -- verificar direto no `cp2k_input.xml`
  gerado por `cp2k.psmp --xml` (ver seção 8, mesma técnica que já pegou o
  bug do `&LIBXC`), não confiar em documentação de outra versão. Por
  segurança, o parser trata isso em Python, não como keyword do `.inp`.
- Confirmar nome do ambiente conda com `pymatgen`/`mp-api` pro
  `job_dos_bands.sh` (`conda activate cp2k_env` está como placeholder).
- Confirmar traceback do `mp-api` (`python3 -c "from mp_api.client import MPRester"`)
  -- suspeita de conflito de versão com `pymatgen` pinado no `vasp_env`; recomendado
  criar um `cp2k_env` conda separado.
- Confirmar e-mail correto em `--mail-user` -- divergência entre
  `roberto.ramos@unesp.br` (jobs de k-mesh/DOS-bandas) e
  `robeerto.aguiar.ramos@gmail.com` (jobs de cutoff); decidir se deve ser o
  mesmo endereço em todos.
- A API key do Materials Project está com default hardcoded em
  `input_generator.py`/`kmesh_stress_convergence.py` (`--api-key`) -- mover para
  exigir `MP_API_KEY` no ambiente (sem default embutido no código) antes de
  compartilhar/versionar estes scripts.

---

## 10. Slabs / energia de superfície (`slab_generator.py`, 2026-09-24)

Depois do Etot (passo 4), cada slab já cortado (ex. POSCAR legacy do VASP) vira um `GEO_OPT`
de célula fixa, consistente com o bulk do CP2K:

```bash
python ~/work_cp2k/scripts/slab_generator.py --slab POSCAR_legacy \
    --bulk-ref <bulk de onde o slab foi cortado> --bulk-cp2k <CIF do Etot CP2K> \
    --cutoff <C> --k-density <n_i*|a_i| da malha de bulk> \
    --project-name <nome> --output-dir <pasta do slab>
sbatch ~/work_cp2k/scripts/job_etot.sh <pasta do slab>
```

- Identifica o plano (hkl) casando os vetores do plano com [uvw] do bulk; recusa fios
  (vácuo em dois eixos) a menos que `--allow-wire`.
- Reescala o plano para a rede do CP2K (mesmos [uvw]) e o eixo normal por d_hkl(CP2K)/d_hkl(ref),
  o que remove a deformação artificial VASP→CP2K.
- Reorienta: a ∥ x, b no plano xy, c ∥ z (normal) com espessura + `--vacuum` (15 Å).
  `SURFACE_DIPOLE_CORRECTION`/`SURF_DIP_DIR Z` (equivalente a `LDIPOL`/`IDIPOL=3`) está ligado
  por padrão e só faz sentido com a normal ao longo de z.
- k-mesh `ceil(k_density/|a|) x ceil(k_density/|b|) x 1`; `ALPHA` de mistura 0.2 (slabs com
  vácuo sofrem mais com charge sloshing).
- Escreve `slab_meta.json` (hkl, n de fórmulas, área, esteq.) para o cálculo de
  gamma = (E_slab − n·E_bulk/fu)/2A. Slab não estequiométrico → gamma depende de mu_O.
