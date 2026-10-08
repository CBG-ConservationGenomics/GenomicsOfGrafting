# GraftingGeneFamilies

Nextflow pipeline for plant gene-family analysis with **OrthoFinder**, functional annotation, and protein percent-identity tables for user-supplied genes of interest.

## Species (default)

| ID  | Species                 | Source                                      |
|-----|-------------------------|---------------------------------------------|
| Sly | *Solanum lycopersicum*  | ITAG4.0 (Sol Genomics) — keeps `Solyc` IDs  |
| Nbe | *Nicotiana benthamiana* | Niben1.0.1 proteins (URL; override if needed)|
| Ath | *Arabidopsis thaliana*  | NCBI RefSeq `GCF_000001735.4`               |
| Cpe | *Cucurbita pepo*        | NCBI RefSeq `GCF_002766695.1`               |

Edit `conf/species.yaml` to change assemblies, switch a species to `source: local`, or add taxa.

## Pipeline stages

1. **Download** proteomes **and CDS** (NCBI Datasets or URL / local path)
2. **QC + rename** proteins (longest isoform; haplotig filter; `Sly|Genus_species|…` IDs)
3. **QC CDS** — keep sequences matching retained proteins; trim transcript ORFs to the protein
4. **BUSCO gate** — protein-mode; species below `--busco_min_complete` (default **90%**) excluded
5. **OrthoFinder** (`-M msa`) — orthogroups + protein MSAs
6. **Select** orthogroups containing genes of interest
7. **Codon alignments** — reverse-translate / project CDS onto OrthoFinder protein MSAs → `05b_codon_alignments/`
8. **Annotate** selected OG proteins (eggNOG, InterProScan, SignalP, DeepTMHMM)
9. **% identity** + polished tables; optional graft phenotype / HyPhy on codon MSAs

AlphaFold is **three-stage**: CPU handoff (`08_alphafold_inputs/`) → GPU fold (`workflows/alphafold.nf`) → CPU structure comparison (`workflows/compare_structures.nf`). Details: [`docs/alphafold_plan.md`](docs/alphafold_plan.md).

Graft phenotype association (outcome matrix → lineage traits, ASR, divergence tests, optional HyPhy/PCOC): [`docs/graft_phenotype.md`](docs/graft_phenotype.md) and `workflows/graft_phenotype.nf`.

## Genes of interest (mixed species)

One entry per line; tomato, Arabidopsis, *N. benthamiana*, and *Cucurbita* IDs can be mixed in the same file:

```text
Solyc01g000570
Ath:AT1G01010
Nbe:Niben101Scf01786g00001
Cpe:XP_023545678.1
```

Formats:

| Format | Example | Notes |
|--------|---------|--------|
| Bare ID | `Solyc01g000570`, `AT1G01010` | Species inferred from ID pattern when possible |
| `Species:gene` | `Ath:AT1G01010` | Recommended; codes: `Sly`, `Ath`, `Nbe`, `Cpe` |
| Tab-separated | `AT1G01010\tAth` | Same as explicit species |

Matching tolerates transcript/isoform suffixes. For ambiguous NCBI accessions (`XP_` / `NP_`), always use an explicit species prefix. Percent identity uses the matched GOI member as the reference by default (`--reference_species_id auto`).

## Run on a remote high-memory server

```bash
# once: install Nextflow + (recommended) Mamba/Conda
conda env create -f conf/environment.yml   # or let Nextflow create it via -profile conda

# set server-side DB / tool paths
export NXF_OPTS='-Xms1g -Xmx8g'

nextflow run main.nf \
  -profile hpc,conda \
  --genes_of_interest data/genes_of_interest.txt \
  --outdir results \
  --eggnog_data_dir /path/to/eggnog_db \
  --interproscan_dir /path/to/interproscan \
  --signalp_bin /path/to/signalp \
  --orthofinder_threads 32
```

### Useful toggles

```bash
--run_eggnog false
--run_interproscan false
--run_signalp false
--run_deeptmhmm false
--reference_species_id Sly

# BUSCO gate
--busco_min_complete 90
--busco_lineage embryophyta_odb10   # or eudicots_odb10
--busco_download_path /path/to/busco_downloads
--busco_offline true               # if lineage DBs already downloaded
--run_busco false                  # skip gate (not recommended)
```

### SLURM

```bash
nextflow run main.nf -profile slurm,conda ...
```

Adjust `queue` / `clusterOptions` in `nextflow.config`.

## Main outputs

| Path | Contents |
|------|----------|
| `results/02_clean_proteomes/*.clean.faa` | QC’d, renamed proteomes |
| `results/02_clean_proteomes/*.id_map.tsv` | original ↔ new ID map |
| `results/02b_busco/busco_summary.tsv` | Per-species BUSCO Complete % + pass/fail |
| `results/02b_busco/gate/*.busco_gate.log` | Gate decision log per species |
| `results/03_orthofinder/` | Full OrthoFinder run |
| `results/05_selected_orthogroups/` | GOI orthogroups + member FASTA |
| `results/06_percent_identity/` | Pairwise + summary %ID |
| `results/07_tables/gene_families_of_interest.tsv` | Polished family table |
| `results/07_tables/gene_families_pctid_by_species.tsv` | Wide %ID by species |
| `results/08_alphafold_inputs/alphafold_inputs/` | Handoff package for the GPU AlphaFold workflow |
| `results/09_structure_comparison/` | RMSD, mapped substitutions, superimposed PDBs |

## Server prerequisites

| Tool | Notes |
|------|--------|
| Nextflow ≥ 23.04 | `curl -s https://get.nextflow.io \| bash` |
| OrthoFinder + DIAMOND | via `conf/environment.yml` |
| eggNOG-mapper + DB | download DB once; pass `--eggnog_data_dir` |
| InterProScan | large install; pass `--interproscan_dir` |
| SignalP 6 | academic license; set `--signalp_bin` |
| DeepTMHMM | often via BioLib; module is a thin wrapper |
| NCBI Datasets CLI | for `source: ncbi` species |

## Local proteomes

If a download URL is flaky (common for *N. benthamiana*), stage files on the server and set:

```yaml
- id: Nbe
  genus: Nicotiana
  species: benthamiana
  source: local
  fasta_path: /data/proteomes/Niben.proteins.fasta
  gff_path: null
```

## Project layout

```text
main.nf
nextflow.config
conf/                 # species.yaml, conda env, HPC resources
modules/local/        # Nextflow process modules
bin/                  # Python helpers (auto on PATH)
data/genes_of_interest.txt
docs/alphafold_plan.md
```
