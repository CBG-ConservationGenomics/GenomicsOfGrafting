# Graft phenotype → sites associated with compatibility

## Goal

Link a **bidirectional scion×rootstock graft-outcome matrix** to sequence change in gene families of interest, while keeping the **OrthoFinder species tree** as the phylogenetic backbone.

## Input matrix

The table is `data/graft_outcomes_11-09-2026.tsv`. A short copy of the same columns is `data/graft_outcomes.example.tsv`.

```text
Rootstock,Scion,success
Solanum lycopersicum,Nicotiana benthamiana,TRUE
Arabidopsis thaliana,Solanum lycopersicum,FALSE
```

`success` is `TRUE` or `FALSE`. Names are scientific names (`Genus species`). A hybrid such as `Solanum lycopersicum × S. pimpinellifolium` is not treated as `Solanum lycopersicum`. The pipeline matches a name to a species id in `conf/species.yaml` only when the binomial is exact.

A run does not have to include every species in that file. `10_graft_phenotype/traits/analyzed_species.tsv` lists the species on the OrthoFinder tree and whether each one has a row in the graft table. Pairs used for divergence are the ones where both species were in the run.

## Analyses implemented

| Step | What it does | Key outputs |
|------|----------------|-------------|
| **Code the phenotype** | Derive pair outcomes, lineage traits (e.g. broad scion), phenotype-shift branches, species distances | `10_graft_phenotype/traits/` |
| **Pairwise divergence** | Scion–rootstock protein %ID per selected OG | `divergence/graft_pair_divergence.tsv` |
| **Divergence vs phenotype** | Mann–Whitney + partial Mantel (protein divergence \| species distance) | `divergence/divergence_vs_phenotype.tsv` |
| **Ancestral reconstruction** | IQ-TREE `-asr`; AA changes on branches; flag phenotype-shift branches | `asr/<OG>/*.asr_substitutions.tsv` |
| **Convergence** | PCOC on the protein alignment, CSUBST on the codon alignment. Images are pulled when the programs are not on PATH | `convergence/<OG>/` |
| **Selection** | HyPhy BUSTED, aBSREL, MEME, RELAX on pipeline codon MSAs. `-profile conda` installs HyPhy from `conf/hyphy.yml` | `selection/<OG>/` |
| **Rank candidates** | Priority score combining the above | `graft_candidate_ranking.tsv` |

In `main.nf`, this ranking is produced **before** the AlphaFold handoff. High-scoring OGs get lower `fold_rank` in `08_alphafold_inputs/` (fold first). Trim the GPU queue with `--af_max_ogs` / `--af_min_priority_score`.

### Interpretation notes (as designed)

- **ASR on shift branches** — strongest *positional* candidates when a substitution arises where the compatibility trait changes.
- **Convergence (PCOC/CSUBST)** — for scattered compatible lineages; phylogeny-aware.
- **Selection (HyPhy)** — supports functional change; does **not** prove a grafting link.
- **Partial Mantel** — pursue a family only if protein divergence predicts failure **beyond** species relatedness (`candidate_beyond_relatedness=yes`).

## How to run

### After `main.nf` (recommended)

```bash
nextflow run workflows/graft_phenotype.nf -profile hpc,conda \
  --graft_matrix data/graft_outcomes_11-09-2026.tsv \
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

HyPhy is installed by `-profile conda` from `conf/hyphy.yml` the first time the selection step runs. PCOC and CSUBST are not in that environment. The convergence step runs a copy already on `PATH`, or pulls `carinerey/pcoc` and `quay.io/biocontainers/csubst:1.14.4--py311h04f3980_0` with Docker or Singularity/Apptainer. The compute node needs network access for that first pull. CSUBST reads the codon alignments in `05b_codon_alignments/`; PCOC reads the protein alignment. Each orthogroup still gets a row in `convergence/<OG>/*.convergence_summary.tsv` when an image cannot be pulled.

### From `main.nf`

```bash
--run_graft_phenotype true \
--graft_matrix data/graft_outcomes_11-09-2026.tsv
```

## Lineage traits

Controlled by:

- `--graft_broad_threshold` (default 0.5) — success rate to call `broad_scion` / `broad_rootstock`
- `--graft_foreground_trait` — `broad_scion` | `broad_rootstock` | `broad_either` (foreground for convergence/selection)

Self-grafts are ignored when computing rates.
