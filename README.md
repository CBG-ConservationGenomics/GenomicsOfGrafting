# GraftingGeneFamilies

Nextflow pipeline for **plant gene-family analysis** in a grafting / scion–rootstock context.

It downloads proteomes and CDS for a set of species, QC’s them, gates annotation quality with BUSCO, builds orthogroups with OrthoFinder, pulls out families containing genes of interest, annotates those proteins, computes protein % identity, prepares AlphaFold inputs, and (optionally) links a bidirectional graft-outcome matrix to candidate sequence changes via ancestral reconstruction, selection, convergence, and divergence tests.

**Status:** v0.1 scaffold — core workflow is implemented and partially run-tested on the HPC (`curie`). Several server-side tool paths and a full end-to-end production run still need to be completed.

---

## Goal (scientific)

Given:

1. A list of plant species  
2. Genes of interest (any of those species’ locus IDs)  
3. Optionally, a scion × rootstock graft-outcome matrix  

Produce:

- Clean, comparable proteomes and orthogroups  
- Gene families of interest with annotations and % identity  
- Codon alignments for selection tests  
- Ranked candidate sites/families associated with graft compatibility  
- A handoff package for AlphaFold structure comparison on a GPU server  

---

## Where we are

### Implemented

| Area | What’s in place |
|------|-----------------|
| **Engine** | Nextflow DSL2 (`main.nf`), conda env, `hpc` / `slurm` / `gpu` profiles |
| **Species config** | `conf/species.yaml` — Sly, Nbe, Ath, Cpe (proteins + CDS URLs/accessions) |
| **Download** | NCBI Datasets or URL/local; auto-decompress gzip; reject HTML/empty files |
| **Protein QC** | Longest isoform; haplotig filters; rename to `Sly\|Genus_species\|origID` + id map |
| **CDS QC** | Match retained proteins; trim transcript ORFs to protein sequence |
| **BUSCO gate** | Protein mode (`embryophyta_odb10`); exclude species below Complete % threshold (default 90%) |
| **OrthoFinder** | MSA mode; emits OGs, alignments, species tree, gene trees |
| **GOI selection** | Mixed-species ID list (`Solyc…`, `Ath:AT…`, `Nbe:…`, `Cpe:XP_…`) |
| **Codon MSAs** | Protein-alignment-guided projection of CDS (pal2nal-style) → `05b_codon_alignments/` |
| **Annotation** | eggNOG-mapper, InterProScan, SignalP, DeepTMHMM on selected OG proteins |
| **% identity** | From OrthoFinder MSAs; reference = matched GOI (`auto`) or fixed species |
| **Report tables** | Polished family table + wide %ID-by-species |
| **AlphaFold** | CPU handoff package; separate GPU ColabFold workflow; CPU structure comparison |
| **Graft phenotype** | Matrix → lineage/pair traits; IQ-TREE ASR; pair divergence + partial Mantel; HyPhy/PCOC/CSUBST hooks; candidate ranking |

### Run / debug notes so far (HPC)

- Fixed Nextflow parse error: `emit: package` → `emit: handoff` (Groovy reserved word).  
- Fixed N. benthamiana 404: use Niben101 proteins under Sol Genomics `annotation/Niben101/`.  
- Fixed gzip-as-FASTA bug: decompress by magic bytes, not by staging filename.  
- Use **`-profile`** (one dash), not `--profile`.  
- Prefer **`-resume`** after module fixes so only invalidated steps re-run.

### Still to do / harden

- Full production run with real GOI list + real graft matrix  
- Confirm server paths for eggNOG DB, InterProScan, SignalP, BUSCO download cache  
- DeepTMHMM / BioLib availability on the server  
- Install PCOC / CSUBST if convergence tests are required  
- Supply or verify CDS completeness for all four species after QC reports  
- AlphaFold: install LocalColabFold on the GPU box; tune `conf/gpu.config`  
- Optional: automate CDS for more species; ChimeraX/PyMOL scripts for AF diffs  
- Example graft matrix is illustrative — replace with experimental outcomes  

---

## Architecture overview

```text
┌──────────────────────────────── main.nf (CPU / HPC) ─────────────────────────────────┐
│                                                                                        │
│  species.yaml ──► DOWNLOAD (protein + CDS + GFF)                                       │
│                      │                                                                 │
│                      ├─► QC proteins ──► BUSCO gate ──► pass proteomes                 │
│                      └─► QC CDS (match protein IDs / trim ORFs)                        │
│                                      │                                                 │
│                                      ▼                                                 │
│                               OrthoFinder (-M msa)                                     │
│                                      │                                                 │
│                    ┌─────────────────┼─────────────────┐                               │
│                    ▼                 ▼                 ▼                               │
│              Select GOI OGs    Species tree      Protein MSAs                          │
│                    │                                   │                               │
│                    ├─► Codon MSAs (CDS ⊕ protein MSA)                                  │
│                    ├─► Annotate (eggNOG / IPS / SignalP / DeepTMHMM)                   │
│                    ├─► % identity + tables                                             │
│                    └─► AlphaFold handoff (08_alphafold_inputs/)                        │
│                                                                                        │
│  optional: --run_graft_phenotype true                                                  │
│            graft matrix + species tree + OGs + codon MSAs                              │
│            → traits, ASR, divergence tests, HyPhy, ranking                             │
└────────────────────────────────────────────────────────────────────────────────────────┘
                    │                                      │
                    ▼                                      ▼
     workflows/alphafold.nf (GPU)           workflows/compare_structures.nf (CPU)
     ColabFold on handoff sequences         RMSD, substitutions, superimposed PDBs
```

Secondary entrypoints:

| Workflow | Role |
|----------|------|
| `main.nf` | Full CPU pipeline |
| `workflows/graft_phenotype.nf` | Phenotype association only (after OrthoFinder exists) |
| `workflows/alphafold.nf` | GPU folding |
| `workflows/compare_structures.nf` | Post-AF structure comparison |

---

## Default species

Configured in [`conf/species.yaml`](conf/species.yaml):

| ID | Species | Proteins | CDS |
|----|---------|----------|-----|
| **Sly** | *Solanum lycopersicum* | ITAG4.0 (Sol Genomics) | ITAG4.0_CDS.fasta |
| **Nbe** | *Nicotiana benthamiana* | Niben101 proteins (gzip) | Niben101 transcripts (ORF-trimmed in QC) |
| **Ath** | *Arabidopsis thaliana* | NCBI `GCF_000001735.4` | NCBI cds_from_genomic |
| **Cpe** | *Cucurbita pepo* | NCBI `GCF_002766695.1` | NCBI cds_from_genomic |

IDs are renamed to:

```text
Sly|Solanum_lycopersicum|Solyc01g000570.3.1
```

so OrthoFinder and downstream tools keep species identity in the sequence header. Original IDs are retained in `*.id_map.tsv`.

To use local files on the server:

```yaml
- id: Nbe
  genus: Nicotiana
  species: benthamiana
  source: local
  fasta_path: /data/proteomes/Niben.proteins.fasta
  cds_path: /data/proteomes/Niben.transcripts.fasta
  gff_path: null
```

---

## Pipeline stages (detail)

### 1–2. Download and protein QC

- Sources: `ncbi` | `url` | `local`  
- Gzip/bzip2 detected by magic bytes and decompressed  
- Keep longest isoform; optionally drop alt-haplotig proteins (GFF/ID patterns)  
- Output: `results/01_raw_proteomes/`, `results/02_clean_proteomes/`

### 3. CDS QC

- Retain CDS matching cleaned proteins via id map  
- If the download is cDNA/transcript (Nben), find the ORF that translates to the protein  
- Output: `results/02_clean_cds/`

### 4. BUSCO gate

- Mode: proteins; default lineage `embryophyta_odb10` (overridable per species)  
- Metric: **Complete %** = Single + Duplicated  
- Default threshold: **90%** (`--busco_min_complete`)  
- Failing species are excluded from OrthoFinder onward  
- Pipeline errors if fewer than `--busco_min_species` (default 2) pass  
- Output: `results/02b_busco/`

### 5. OrthoFinder

- `-M msa` so MSAs exist for %ID and codon projection  
- Emits orthogroups, MSAs, rooted species tree, gene trees  
- Output: `results/03_orthofinder/`

### 6. Genes of interest → orthogroups

Input file (`data/genes_of_interest.txt` or your list), one entry per line:

```text
Solyc01g000570
Ath:AT1G01010
Nbe:Niben101Scf01786g00001
Cpe:XP_023545678.1
```

| Format | Example |
|--------|---------|
| Bare ID | `Solyc01g000570` (species inferred when possible) |
| `Species:gene` | `Ath:AT1G01010` — preferred for ambiguous IDs |
| Tab-separated | `AT1G01010\tAth` |

Output: `results/05_selected_orthogroups/`

### 7. Codon alignments

For each selected OG, project CDS onto the OrthoFinder protein MSA:

```text
protein MSA:  M K T - F V
CDS:          ATGAAA ACT   TTCGTT
codon MSA:    ATGAAA ACT --- TTCGTT
```

Output: `results/05b_codon_alignments/OG*.codon.fna`  
Used by HyPhy when `--run_hyphy true`.

### 8. Functional annotation

Run on **selected OG proteins only** (not whole proteomes):

- eggNOG-mapper  
- InterProScan  
- SignalP  
- DeepTMHMM  

Toggle with `--run_eggnog false`, etc. Point DB/tool installs with `--eggnog_data_dir`, `--interproscan_dir`, `--signalp_bin`.

### 9. % identity and tables

- Pairwise and summary %ID vs matched GOI (or `--reference_species_id Sly`)  
- Main table: `results/07_tables/gene_families_of_interest.tsv`  
- Wide table: `gene_families_pctid_by_species.tsv`

### 10. AlphaFold (optional, multi-machine)

1. **CPU** (`main.nf`): writes `results/08_alphafold_inputs/alphafold_inputs/`  
   (`sequences/`, `manifest.tsv`, `pairs.tsv`)  
2. **GPU** (`workflows/alphafold.nf -profile gpu`): ColabFold  
3. **CPU** (`workflows/compare_structures.nf`): CA RMSD, substitution list + pLDDT, superimposed PDBs → `09_structure_comparison/`  

See [`docs/alphafold_plan.md`](docs/alphafold_plan.md).

### 11. Graft phenotype association (optional)

Enable with `--run_graft_phenotype true` or run `workflows/graft_phenotype.nf` after OrthoFinder.

Input matrix (bidirectional): `data/graft_outcomes.example.tsv`

```text
scion	rootstock	outcome
Sly	Nbe	compatible
Nbe	Sly	compatible
Sly	Ath	incompatible
```

Outcomes: `compatible` | `incompatible` | `delayed_failure`

Analyses:

| Analysis | Role |
|----------|------|
| Lineage + pair traits | “Broad scion/rootstock,” reciprocity, phenotype-shift branches |
| IQ-TREE `-asr` | AA changes; flag those on phenotype-shift branches |
| Pair divergence + partial Mantel | Does protein divergence predict failure **beyond** species relatedness? |
| HyPhy (BUSTED, aBSREL, MEME, RELAX) | Selection on codon MSAs (foreground = broad-grafting tips) |
| PCOC / CSUBST | Convergence (if installed) |
| Ranking | `graft_candidate_ranking.tsv` |

Details: [`docs/graft_phenotype.md`](docs/graft_phenotype.md).

---

## How to run

### Main pipeline (HPC)

```bash
conda env create -f conf/environment.yml   # or rely on -profile conda
export NXF_OPTS='-Xms1g -Xmx8g'

nextflow run main.nf \
  -profile hpc,conda \
  --genes_of_interest data/grafting_candidates.txt \
  --outdir results \
  --eggnog_data_dir /home/general/Databases/Annotation/eggnog/ \
  --interproscan_dir /home/general/Databases/Annotation/interproscan \
  --signalp_bin /path/to/signalp \
  --busco_download_path /path/to/busco_downloads \
  --orthofinder_threads 32
```

Resume after a fix:

```bash
nextflow run main.nf -profile hpc,conda ... -resume
```

### With graft phenotype + HyPhy

```bash
nextflow run main.nf -profile hpc,conda \
  ... \
  --run_graft_phenotype true \
  --graft_matrix data/graft_outcomes.tsv \
  --build_codon_alignments true \
  --run_hyphy true
```

### Standalone graft phenotype (after OrthoFinder exists)

```bash
nextflow run workflows/graft_phenotype.nf -profile hpc,conda \
  --graft_matrix data/graft_outcomes.tsv \
  --species_tree results/03_orthofinder/OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt \
  --og_list results/05_selected_orthogroups/selected_orthogroups.txt \
  --alignments_dir results/03_orthofinder/OrthoFinder/Results_*/MultipleSequenceAlignments \
  --outdir results
```

### AlphaFold on the GPU server

```bash
rsync -av results/08_alphafold_inputs/alphafold_inputs/ gpu:/data/handoff/alphafold_inputs/

nextflow run workflows/alphafold.nf -profile gpu \
  --alphafold_inputs /data/handoff/alphafold_inputs \
  --outdir alphafold_results
```

Then compare structures on CPU:

```bash
nextflow run workflows/compare_structures.nf -profile hpc,conda \
  --alphafold_inputs results/08_alphafold_inputs/alphafold_inputs \
  --alphafold_models alphafold_results/models \
  --alignments_dir results/03_orthofinder/OrthoFinder/Results_*/MultipleSequenceAlignments \
  --outdir results
```

---

## Important parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `--genes_of_interest` | `data/genes_of_interest.txt` | GOI list |
| `--species_yaml` | `conf/species.yaml` | Species / download config |
| `--busco_min_complete` | `90` | Min Complete % to keep a species |
| `--busco_lineage` | `embryophyta_odb10` | Default BUSCO lineage |
| `--reference_species_id` | `auto` | %ID reference (matched GOI or `Sly`/`Ath`/…) |
| `--build_codon_alignments` | `true` | Build codon MSAs from CDS + protein MSAs |
| `--run_graft_phenotype` | `false` | Run phenotype association at end of `main.nf` |
| `--graft_matrix` | example TSV | Scion×rootstock outcomes |
| `--run_hyphy` | `false` | Selection tests (needs codon MSAs) |
| `--run_asr` | `true` | IQ-TREE ancestral reconstruction |
| `--run_convergence` | `true` | PCOC/CSUBST if installed |
| `--prepare_alphafold` | `true` | Write AF handoff package |
| `--orthofinder_threads` | `16` | OrthoFinder CPUs |

Annotation / DB paths: `--eggnog_data_dir`, `--interproscan_dir`, `--signalp_bin`, `--busco_download_path`.

---

## Main outputs

| Path | Contents |
|------|----------|
| `01_raw_proteomes/` | Downloaded proteins, CDS, GFF |
| `02_clean_proteomes/` | QC’d proteins + id maps |
| `02_clean_cds/` | QC’d CDS matching proteins |
| `02b_busco/` | BUSCO summaries + pass/fail gate logs |
| `03_orthofinder/` | Orthogroups, MSAs, species/gene trees |
| `05_selected_orthogroups/` | GOI orthogroups + member FASTA |
| `05b_codon_alignments/` | Codon MSAs for selected OGs |
| `04_annotation/` | eggNOG / InterProScan / SignalP / DeepTMHMM |
| `06_percent_identity/` | Pairwise + summary %ID |
| `07_tables/` | Polished family tables |
| `08_alphafold_inputs/` | GPU handoff package |
| `09_structure_comparison/` | RMSD, substitutions, superimposed PDBs |
| `10_graft_phenotype/` | Traits, ASR, divergence tests, ranking |

---

## Project layout

```text
GraftingGeneFamilies/
├── main.nf                          # Primary CPU pipeline
├── nextflow.config                  # Params + profiles
├── conf/
│   ├── species.yaml                 # Species, URLs, accessions, BUSCO lineage
│   ├── environment.yml              # Conda dependencies
│   ├── base.config / hpc.config / gpu.config
├── modules/local/                   # Nextflow process modules
├── subworkflows/local/
│   └── graft_phenotype.nf
├── workflows/
│   ├── alphafold.nf                 # GPU folding
│   ├── compare_structures.nf        # Post-AF comparison
│   └── graft_phenotype.nf           # Standalone phenotype association
├── bin/                             # Python helpers (on PATH in Nextflow)
├── data/
│   ├── genes_of_interest.txt
│   ├── grafting_candidates.txt      # Your GOI list on the server
│   └── graft_outcomes.example.tsv
├── docs/
│   ├── alphafold_plan.md
│   └── graft_phenotype.md
└── assets/empty.tsv
```

---

## Server prerequisites

| Component | Notes |
|-----------|--------|
| Nextflow ≥ 23.04 | Installed on `curie` |
| Conda env | `conf/environment.yml` (OrthoFinder, BUSCO, IQ-TREE, HyPhy, Biopython, …) |
| eggNOG DB | `--eggnog_data_dir` |
| InterProScan | `--interproscan_dir` |
| SignalP 6 | Licensed binary; `--signalp_bin` |
| BUSCO lineages | Prefer shared `--busco_download_path` |
| DeepTMHMM | BioLib or local install |
| LocalColabFold | GPU server only |
| PCOC / CSUBST | Optional; convergence steps skip cleanly if missing |
| USalign / TMalign | Optional TM-score in structure comparison |

---

## Design choices (short)

- **Nextflow** for HPC resume, profiles, and process isolation  
- **ITAG / Niben community IDs** so GOI lists stay human-readable  
- **BUSCO gate** before OrthoFinder so bad annotations don’t pollute families  
- **Annotate selected OGs only** to save InterProScan/eggNOG time  
- **Codon MSAs from protein MSAs** so HyPhy uses the same homology hypothesis as OrthoFinder  
- **AlphaFold split CPU/GPU** so OrthoFinder and folding don’t fight for resources  
- **Phenotype tests control for phylogeny** (shift branches, partial Mantel, HyPhy foregrounds) rather than raw residue counting  

---

## Suggested next steps

1. Replace example GOIs / graft matrix with experimental data  
2. Point all DB/tool params at real paths on `curie`  
3. Complete a `-resume` run through tables + codon alignments  
4. Turn on `--run_graft_phenotype true` once the matrix is real  
5. Sync AF handoff to the GPU server and fold prioritized pairs  
6. Iterate on candidates that are (a) on phenotype-shift branches, (b) significant beyond species distance, and (c) structurally plausible in AF comparisons  

Questions or design changes: adjust `conf/species.yaml` and `nextflow.config` first — most behavior is parameterized there.
