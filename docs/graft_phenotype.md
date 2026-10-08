# Graft phenotype → sites associated with compatibility

## Goal

Link a **bidirectional scion×rootstock graft-outcome matrix** to sequence change in gene families of interest, while keeping the **OrthoFinder species tree** as the phylogenetic backbone.

## Input matrix

Long TSV (preferred):

```text
scion	rootstock	outcome
Sly	Nbe	compatible
Nbe	Sly	compatible
Sly	Ath	incompatible
```

Outcomes: `compatible` | `incompatible` | `delayed_failure`  
Species IDs must match pipeline ids (`Sly`, `Ath`, `Nbe`, `Cpe`).

Example: `data/graft_outcomes.example.tsv`

## Analyses implemented

| Step | What it does | Key outputs |
|------|----------------|-------------|
| **Code the phenotype** | Derive pair outcomes, lineage traits (e.g. broad scion), phenotype-shift branches, species distances | `10_graft_phenotype/traits/` |
| **Pairwise divergence** | Scion–rootstock protein %ID per selected OG | `divergence/graft_pair_divergence.tsv` |
| **Divergence vs phenotype** | Mann–Whitney + partial Mantel (protein divergence \| species distance) | `divergence/divergence_vs_phenotype.tsv` |
| **Ancestral reconstruction** | IQ-TREE `-asr`; AA changes on branches; flag phenotype-shift branches | `asr/<OG>/*.asr_substitutions.tsv` |
| **Convergence** | PCOC / CSUBST if installed (else skipped with log) | `convergence/<OG>/` |
| **Selection** | HyPhy BUSTED, aBSREL, MEME, RELAX on pipeline codon MSAs | `selection/<OG>/` |
| **Rank candidates** | Priority score combining the above | `graft_candidate_ranking.tsv` |

### Interpretation notes (as designed)

- **ASR on shift branches** — strongest *positional* candidates when a substitution arises where the compatibility trait changes.
- **Convergence (PCOC/CSUBST)** — for scattered compatible lineages; phylogeny-aware.
- **Selection (HyPhy)** — supports functional change; does **not** prove a grafting link.
- **Partial Mantel** — pursue a family only if protein divergence predicts failure **beyond** species relatedness (`candidate_beyond_relatedness=yes`).

## How to run

### After `main.nf` (recommended)

```bash
nextflow run workflows/graft_phenotype.nf -profile hpc,conda \
  --graft_matrix data/graft_outcomes.tsv \
  --species_tree results/03_orthofinder/OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt \
  --og_list results/05_selected_orthogroups/selected_orthogroups.txt \
  --alignments_dir results/03_orthofinder/OrthoFinder/Results_*/MultipleSequenceAlignments \
  --outdir results \
  --run_asr true \
  --run_convergence true \
  --run_hyphy false
```

### With HyPhy (codon alignments)

`main.nf` downloads CDS with proteomes, QC/renames them to match protein IDs, then builds
**protein-MSA-guided codon alignments** (pal2nal-style) into `results/05b_codon_alignments/`.

```bash
--build_codon_alignments true \
--run_hyphy true
```

Optional override with external codon MSAs: `--codon_alignments_dir path/to/codon_alns`.

### From `main.nf`

```bash
--run_graft_phenotype true \
--graft_matrix data/graft_outcomes.tsv
```

## Lineage traits

Controlled by:

- `--graft_broad_threshold` (default 0.5) — success rate to call `broad_scion` / `broad_rootstock`
- `--graft_foreground_trait` — `broad_scion` | `broad_rootstock` | `broad_either` (foreground for convergence/selection)

Self-grafts are ignored when computing rates.
