# Pipeline CP2K — TiO2 (reescrito do zero)

Documento de referência. Substitui integralmente o fluxo anterior, que foi
montado portando convenções do VASP e continha erros estruturais.

Sistema alvo: TiO2 célula unitária (≤10 átomos), não-magnético, PBE →
estrutura eletrônica.

---

## PARTE 1 — Achados: o que estava errado e por quê

Todos verificados na documentação oficial do CP2K (manual.cp2k.org,
cp2k.org/exercises) e no paper do pacote (arXiv:2508.15559, arXiv:2003.03868).

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

---

## PARTE 2 — Pipeline novo

Quatro etapas. Cada uma com um gerador Python e um job SLURM.

```
[0] converge CUTOFF + REL_CUTOFF   (método 2D do manual, MAX_SCF 1)
         ↓
[1] converge k-mesh                (CELL_OPT completo por malha)
         ↓
[2] CELL_OPT de produção           (PASSO ÚNICO — célula + íons)
         ↓
    ├── [3a] BANDAS   → k-mesh + &PRINT%BAND_STRUCTURE
    └── [3b] PDOS     → supercélula, Gamma-only, OT, &PRINT%PDOS
```

### Etapa 0 — CUTOFF e REL_CUTOFF (método oficial)

- `RUN_TYPE ENERGY`, `MAX_SCF 1`, `PRINT_LEVEL MEDIUM`, Gamma-only
- Fase A: varrer CUTOFF (200→1000 Ry) com REL_CUTOFF fixo em 60
- Fase B: varrer REL_CUTOFF (10→100 Ry) com o CUTOFF escolhido
- Parser extrai **energia total E a distribuição de Gaussianas por grid**
- Critério: menor CUTOFF que usa o grid mais fino, com a maioria das
  Gaussianas nos grids grossos

### Etapa 1 — k-mesh

- `CELL_OPT` completo por malha (agora barato: passo único, ≤10 átomos)
- Malhas 2×2×2 → 12×12×12
- Critério duplo: ΔE (meV/átomo) **e** ΔVolume (%) vs. a malha mais densa
- Com `DIAGONALIZATION` + `SMEAR` (k≠Gamma não aceita OT)

### Etapa 2 — CELL_OPT de produção (passo único)

- `OPTIMIZER BFGS`, `MAX_FORCE [eV*angstrom^-1] 0.010`
- `CELL_REF` = 1.15× os vetores da célula de entrada
- `KEEP_SPACE_GROUP` opcional (preserva simetria durante a otimização)
- `STRESS_TENSOR ANALYTICAL`
- `EXTRAPOLATION USE_GUESS`, `&SMEAR`, `PARALLEL_GROUP_SIZE -1`
- Escreve `.wfn` explicitamente (`&SCF%PRINT%RESTART ON`)

### Etapa 3a — Estrutura de bandas

- `RUN_TYPE ENERGY` sobre a estrutura relaxada
- Mesmo k-mesh da etapa 2 (bandas NÃO precisam de malha mais densa)
- `ADDED_MOS 10` no `&SCF` e no `&BAND_STRUCTURE`
- Caminho de alta simetria via pymatgen `HighSymmKpath`
- **Atenção:** usar `kpath.prim` (célula padronizada) no `&CELL`/`&COORD`,
  senão as coordenadas do caminho não correspondem à rede simulada

### Etapa 3b — PDOS (rota da supercélula)

- Supercélula 2×2×2 (ou o que a convergência indicar) da célula relaxada
- **Gamma-only** (`&KPOINTS` ausente ou `SCHEME NONE`)
- `&OT` (volta a ser possível, porque é Gamma)
- `&PRINT%PDOS` com `COMPONENTS` e `NLUMO -1`
- Alargamento gaussiano (equivalente ao `NEDOS=6000`) feito **em Python**,
  no parser, sobre os `.pdos` crus — não como keyword do input (ver 1.9)

---

## PARTE 3 — Pendências a verificar no cluster

1. `&DFT%PRINT%DOS` existe na build 2026.1? (teste rápido; se abortar → rota 3b)
2. RI-HFX com k-points funciona na 2026.1? (decide se bandas HSE06 são viáveis)
3. Marcador exato de convergência no output (`OPTIMIZATION COMPLETED`?)
4. Nome do ambiente conda com pymatgen/mp-api
5. Tamanho da supercélula para o PDOS (requer teste de convergência próprio)

---

## PARTE 4 — Referências

- manual.cp2k.org/trunk/methods/dft/cutoff.html — convergência CUTOFF/REL_CUTOFF
- manual.cp2k.org/trunk/methods/optimization/geometry.html — CELL_REF, geo/cell opt
- manual.cp2k.org/trunk/methods/dft/k-points.html — k-points, wavefunctions complexas
- manual.cp2k.org/trunk/methods/electronic_structure/dos.html — DOS/PDOS, limitação de versão
- cp2k.org/exercises:common:bs — exemplo oficial de bandas (WO3)
- cp2k.org/exercises:2017_uzh_cmest:pdos — limitação PDOS + k-points
- github.com/cp2k/cp2k/issues/4854 — status de suporte a k-points
- arXiv:2508.15559 — *The CP2K Program Package Made Simple*
- arXiv:2003.03868 — *CP2K/Quickstep*

---

## PARTE 5 — Adendo 2026-09-22: `&HF/&MEMORY`, `T_C_G_DATA` e a interface nova de PDOS

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
