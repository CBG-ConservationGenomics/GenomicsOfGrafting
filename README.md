q# GraftingGeneFamilies

Nextflow pipeline for **plant gene-family analysis** in a grafting / scion–rootstock context.

It downloads proteomes and CDS for a set of species, QC’s them, gates annotation quality with BUSCO, builds orthogroups with OrthoFinder, pulls out families containing genes of interest, annotates those proteins, and computes protein % identity. Optionally it then links a bidirectional graft-outcome matrix to candidate sequence changes via ancestral reconstruction, selection, convergence, and divergence tests. The AlphaFold handoff is the last CPU step, so the GPU queue can be limited to families that ranking (or % identity, if phenotype is off) actually selects.

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
- A selective AlphaFold handoff, written last, for structure comparison on a GPU server  

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
| **Candidate CDS check** | Optional characterized CDS FASTA; confirms the annotation CDS and protein for each resolved gene |
| **Codon MSAs** | Protein-alignment-guided projection of CDS (pal2nal-style) → `05b_codon_alignments/` |
| **Annotation** | eggNOG-mapper, InterProScan, SignalP, DeepTMHMM on selected OG proteins |
| **% identity** | From OrthoFinder MSAs; reference = matched GOI (`auto`) or fixed species |
| **Report tables** | Polished family table + wide %ID-by-species |
| **Graft phenotype** | Matrix → lineage/pair traits; IQ-TREE ASR; pair divergence + partial Mantel; HyPhy/PCOC/CSUBST hooks; candidate ranking |
| **AlphaFold** | Last CPU step: handoff limited to prioritized families; separate GPU ColabFold workflow; CPU structure comparison |

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
│                    ├─► optional characterized CDS check                                │
│                    ├─► Codon MSAs (CDS ⊕ protein MSA)                                  │
│                    ├─► Annotate (eggNOG / IPS / SignalP / DeepTMHMM)                   │
│                    └─► % identity + tables                                             │
│                                      │                                                 │
│                                      ▼                                                 │
│                    optional graft phenotype (ASR / HyPhy / ranking)                    │
│                                      │                                                 │
│                                      ▼                                                 │
│                    AlphaFold handoff (fold_rank from ranking)                          │
│                                      │                                                 │
│                                      ▼                                                 │
│                    11_curated/index.html  ← browsable report of the tables             │
└────────────────────────────────────────────────────────────────────────────────────────┘
                    │                                      │
                    ▼                                      ▼
     workflows/alphafold.nf (GPU)           workflows/compare_structures.nf (CPU)
     ColabFold (high fold_rank first)       RMSD, substitutions, superimposed PDBs
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
| **Cpe** | *Cucurbita pepo* | NCBI `GCF_002806865.2` | NCBI cds_from_genomic |

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
- Each download writes `<species>.provenance.tsv`: the accession or URL, the file taken from the package, the organism and annotation the source reported, the protein count, and an MD5  
- An NCBI package whose organism does not match the configured genus stops the download  
- Gzip/bzip2 detected by magic bytes and decompressed  
- Keep longest isoform; optionally drop alt-haplotig proteins (GFF/ID patterns)  
- Output: `results/01_raw_proteomes/`, `results/02_clean_proteomes/`
- If `<species>.proteins.faa`, `.cds.fna`, `.gff3`, and `.provenance.tsv` are already in `01_raw_proteomes/`, that species is not downloaded again. `--reuse_downloads false` fetches them again.
- If `<species>.clean.faa` and `<species>.id_map.tsv` are already in `02_clean_proteomes/`, protein QC is skipped. `--reuse_proteome_qc false` runs it again.

### 3. CDS QC

- Retain CDS matching cleaned proteins via id map  
- If the download is cDNA/transcript (Nben), find the ORF that translates to the protein  
- Output: `results/02_clean_cds/`
- If `<species>.cds_qc_report.tsv` is already in `02_clean_cds/`, CDS QC is skipped and `<species>.clean.cds.fna` is reused when it exists. `--reuse_cds_qc false` runs it again.

### 4. BUSCO gate

- Mode: proteins; default lineage `embryophyta_odb10` (overridable per species)  
- Metric: **Complete %** = Single + Duplicated  
- Default threshold: **90%** (`--busco_min_complete`)  
- Failing species are excluded from OrthoFinder onward  
- Pipeline errors if fewer than `--busco_min_species` (default 2) pass  
- Output: `results/02b_busco/`
- If `02b_busco/<species>/<species>.busco_metrics.tsv` and `busco_status.txt` are already present, BUSCO is not run again for that species. The saved pass/fail is reused, so change `--busco_min_complete` only with `--reuse_busco false`.

### 5. OrthoFinder

- `-M msa` so MSAs exist for %ID and codon projection  
- Emits orthogroups, MSAs, rooted species tree, gene trees  
- Output: `results/03_orthofinder/`

### 6. Genes of interest → orthogroups

Input file (`data/genes_of_interest.txt` or your list). Each line is the functionally characterized gene and the species it was described in:

```text
WOX4	Solanum lycopersicum
AT1G46480	Arabidopsis thaliana
XP_023545678.1	Cucurbita pepo
```

| Column | Meaning |
|--------|---------|
| 1 | Gene symbol (`WOX4`) or locus ID (`Solyc04g078650`, `AT1G46480`) |
| 2 | Genus and species (`Solanum lycopersicum`). A short id such as `Sly` is also accepted |

Symbols and locus IDs are normalized before they are compared with the annotation. `Solyc04g078650` matches the protein `Solyc04g078650.4.1`. A dotted BAC gene name such as `F8D20.60` stays whole, so the clone name is not shared by every gene on that BAC. An AGI id such as `AT4G35550` stays on the gene whose `locus_tag` is that id. A species prefix is removed (`SlWOX4` is searched as `WOX4`, and `GH9B3` matches the synonym `ATGH9B3`). A symbol matches a GFF gene symbol when the annotation has one. It also matches a functional description that uses the expanded family name and the same number: ITAG stores tomato WOX4 as “wuschel related homobox 4”, with no `WOX4` token. A `-like` name is not matched to a different number (`WRKY70-like` does not become WRKY80). The protein that matched is stored as `reference_member` and is the reference for % identity, the family table, and AlphaFold pairs. How each name resolved is written to `gene_id_crosswalk.tsv`. If several literature genes fall in one orthogroup, the first one in this file is the reference. A characterized gene that OrthoFinder left unassigned is still selected, under an `UNASSIGNED_` id, so annotation can run on that protein.

#### Optional characterized CDS check

Gene names can resolve to the wrong locus. To check that, pass a FASTA of the CDS that was characterized, one sequence per candidate, from the species it was described in:

```text
>WOX13	Arabidopsis thaliana AT4G35550.1 TAIR10
ATGG...
```

The gene name and the rest of the header are separated by a tab. After the tab, spaces separate the genus and species, the locus id, and the annotation name. The species is `Arabidopsis thaliana`. `AT4G35550.1` and `TAIR10` are not part of the species. The record matches a candidate listed as `WOX13` or as `AT4G35550`. A species is required when the same symbol is listed for more than one species. `data/candidate_cds.example.fna` shows this layout; its nucleotides are placeholders, not real genes.

The pipeline translates that CDS and aligns it to the cleaned annotation protein chosen for the gene. It also aligns the nucleotides to the cleaned annotation CDS. NCBI `cds_from_genomic` records (`lcl|NC_..._cds_NP_195280.1_1`, with `[protein_id=NP_195280.1]`) are renamed to the same cleaned id as that protein. A terminal stop codon is ignored. Both identities default to a 98% minimum (`--candidate_cds_min_protein_identity`, `--candidate_cds_min_cds_identity`). Lower the CDS threshold when the published cDNA and the genome annotation differ by synonymous substitutions.

The report is `results/05_selected_orthogroups/candidate_cds_check.tsv`.

| status | Meaning |
|--------|---------|
| `match` | Annotation CDS and protein both meet the thresholds |
| `cds_mismatch` | Protein matches; CDS identity is below the threshold |
| `mismatch` | Protein identity is below the threshold. This is still a mismatch when the annotation CDS is missing |
| `missing_annotation_cds` / `missing_annotation_protein` | The resolved id is absent from the cleaned sequences |
| `no_reference` | This candidate has no sequence in the FASTA |
| `unmatched_header` / `ambiguous_header` / `duplicate_header` | The FASTA record could not be tied to one candidate |
| `bad_cds` | The characterized sequence is not a translatable CDS |
| `unresolved` | The gene name was not found in the annotation |

With `--candidate_cds_strict true` (the default), any of those except a bare `unresolved` name stops the run before functional annotation. `no_reference` also stops the run when `--candidate_cds_require_all true` (the default). Set `--candidate_cds_strict false` to write the report and continue.

Output: `results/05_selected_orthogroups/`

### 7. Codon alignments

For each selected OG, project CDS onto the OrthoFinder protein MSA:

```text
protein MSA:  M K T - F V
CDS:          ATGAAA ACT   TTCGTT
codon MSA:    ATGAAA ACT --- TTCGTT
```

Output: `results/05b_codon_alignments/OG*.codon.fna` and `codon_alignment_report.tsv`.  
Alignment headers are matched to cleaned CDS even when OrthoFinder has rewritten `|` as `_`, kept only the accession, or stored an internal id (`0_42`) in `SequenceIDs.txt`. An `UNASSIGNED_` orthogroup has no protein MSA, so its report status is `unassigned_no_alignment`.  
Used by HyPhy when `--run_hyphy true`.

### 8. Functional annotation

Run on **selected OG proteins only** (not whole proteomes):

- eggNOG-mapper  
- InterProScan  
- SignalP  
- DeepTMHMM  

Toggle with `--run_eggnog false`, etc. Point DB/tool installs with `--eggnog_data_dir`, `--interproscan_dir`, `--signalp_bin`.

### 9. % identity and tables

- Pairwise identity of every protein in the alignment, not only the species mean. Square matrices are in `results/06_percent_identity/identity_matrices/`. The characterized protein is the first row and column.  
- Summary %ID versus that protein is still in `summary_percent_identity.tsv`  
- Main table: `results/07_tables/gene_families_of_interest.tsv`  
- Wide table: `gene_families_pctid_by_species.tsv`. `best_pct_id_*` is the closest protein from that species to the reference. `mean_pct_id_*` is still the average of every homolog from that species.

Stages 10 and 11 always run in this order on the CPU: phenotype ranking first, AlphaFold handoff last. Folding itself stays on the GPU.

### 10. Graft phenotype association (optional; before the AlphaFold handoff)

Enable with `--run_graft_phenotype true` or run `workflows/graft_phenotype.nf` after OrthoFinder.

Runs **before** the AlphaFold handoff so ASR / divergence / HyPhy ranking can prioritize what to fold.

Input matrix: `data/graft_outcomes_11-09-2026.tsv`

```text
Rootstock,Scion,success
Solanum lycopersicum,Nicotiana benthamiana,TRUE
Arabidopsis thaliana,Solanum lycopersicum,FALSE
```

`success` is `TRUE` or `FALSE`. Names are scientific names. Only the species in the current run are compared, and `10_graft_phenotype/traits/analyzed_species.tsv` records each of those species and whether the graft table has an entry for it.

Analyses:

| Analysis | Role |
|----------|------|
| Lineage + pair traits | “Broad scion/rootstock,” reciprocity, phenotype-shift branches |
| IQ-TREE `-asr` | AA changes; flag those on phenotype-shift branches |
| Pair divergence + partial Mantel | Does protein divergence predict failure **beyond** species relatedness? |
| HyPhy (BUSTED, aBSREL, MEME, RELAX) | Selection on codon MSAs (foreground = broad-grafting tips). `-profile conda` installs it from `conf/hyphy.yml` |
| PCOC / CSUBST | Convergence. The step pulls `carinerey/pcoc` and the CSUBST Biocontainer when those programs are not already on PATH |
| Ranking | `graft_candidate_ranking.tsv` → feeds AlphaFold handoff priority |

Details: [`docs/graft_phenotype.md`](docs/graft_phenotype.md).

### 11. AlphaFold handoff (optional; last analysis step)

1. **CPU** (`main.nf`, after graft phenotype when that flag is on): writes `results/08_alphafold_inputs/alphafold_inputs/`  
   (`sequences/`, `manifest.tsv` sorted by `fold_rank`, `og_priority.tsv`, `pairs.tsv`). The `08_` folder name is unchanged; it is published after `10_graft_phenotype/`.  
2. **GPU** (`workflows/alphafold.nf -profile gpu`): ColabFold (fold high `fold_rank` first)  
3. **CPU** (`workflows/compare_structures.nf`): CA RMSD, substitution list + pLDDT → `09_structure_comparison/`  

Shrink the GPU queue with `--af_max_ogs 20` or `--af_min_priority_score 5` (requires `--run_graft_phenotype true`).

See [`docs/alphafold_plan.md`](docs/alphafold_plan.md).

**Functional domains on structures:** map SignalP / InterProScan / DeepTMHMM / N-glycosylation sequons onto AF residue numbers; flag scion–rootstock substitutions inside those features; prepare AF-Multimer jobs for scion ligand × rootstock receptor. Also: multimer interface PAE/contacts, defense surface RSA, curated sites via OrthoFinder MSA. See [`docs/structure_feature_mapping.md`](docs/structure_feature_mapping.md).

### 12. Curated report

Written at the end of `main.nf` to `results/11_curated/`. Open `index.html` in a browser to search, sort, and expand each orthogroup. The OrthoFinder species tree is drawn at the top of the family table. Open a family to see its pairwise identity heatmap beside that family's OrthoFinder gene tree. The counts at the top are literature genes placed in an orthogroup, species included in the job, species whose proteins were downloaded from NCBI or a URL, and the time from job start until the report was written. `genomes.html` shows proteins and CDS per genome, a BUSCO composition bar, and a BUSCO heatmap. `protein_alignments.html` shows the OrthoFinder amino-acid alignment for each selected orthogroup. Checkboxes show or hide InterPro domains, SignalP signal peptides, DeepTMHMM transmembrane helices, and N-glycosylation sequons when those predictions exist. `codon_alignments.html` shows each CDS alignment projected onto that protein MSA, with the same checkboxes. `structures.html` opens AlphaFold models in Mol*. Choose the model files, or the superimposed PDBs plus the reference model, and color them by pLDDT, by those same domain tracks, or by amino-acid differences. The `tables/` folder holds copies of the family table, the wide % identity table, and, when those steps ran, the graft ranking, BUSCO summary, codon report, and AlphaFold queue.

---

## How to run

### Server paths and options

Copy `conf/server.config.example` to `conf/server.config` on the server. That example already has the eggNOG, InterProScan, and SignalP paths, 40 OrthoFinder threads, the candidate CDS file, graft phenotype, and HyPhy. The pipeline reads `conf/server.config` when the file exists. It is not committed. A flag on the command line still overrides the file, which is how each run picks its own `--outdir`.

```bash
cp conf/server.config.example conf/server.config
```

### Main pipeline (HPC)

```bash
conda env create -f conf/environment.yml   # or rely on -profile conda
export NXF_OPTS='-Xms1g -Xmx8g'

screen -L nextflow run main.nf -profile hpc,conda --outdir test_out8
```

`conf/server.config` supplies the gene list, `data/candidate_cds.example.fna`, the graft table, HyPhy, and the database paths. Resume the same folder with `-resume`:

```bash
screen -L nextflow run main.nf -profile hpc,conda --outdir test_out8 -resume
```

### Standalone graft phenotype (after OrthoFinder exists)

```bash
nextflow run workflows/graft_phenotype.nf -profile hpc,conda \
  --graft_matrix data/graft_outcomes_11-09-2026.tsv \
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
| `--genes_of_interest` | `data/genes_of_interest.txt` | `gene<TAB>Genus species` for each characterized gene |
| `--candidate_cds` | unset | Characterized CDS FASTA used to confirm each resolved annotation sequence |
| `--candidate_cds_min_protein_identity` | `98` | Minimum identity between the CDS translation and the annotation protein |
| `--candidate_cds_min_cds_identity` | `98` | Minimum identity between the characterized CDS and the annotation CDS |
| `--candidate_cds_require_all` | `true` | Every resolved candidate must have a reference CDS |
| `--candidate_cds_strict` | `true` | Stop the run when the check fails |
| `--species_yaml` | `conf/species.yaml` | Species / download config |
| `--reuse_downloads` | `true` | Skip download when that species is already in `01_raw_proteomes/` |
| `--reuse_proteome_qc` | `true` | Skip protein QC when `02_clean_proteomes/` already has the cleaned FASTA and id map |
| `--reuse_cds_qc` | `true` | Skip CDS QC when `02_clean_cds/` already has the QC report |
| `--reuse_busco` | `true` | Skip BUSCO when that species already has metrics and a pass/fail status |
| `--busco_min_complete` | `90` | Min Complete % to keep a species |
| `--busco_lineage` | `embryophyta_odb10` | Default BUSCO lineage |
| `--reference_species_id` | `auto` | Fallback species if the characterized protein is absent from the alignment |
| `--build_codon_alignments` | `true` | Build codon MSAs from CDS + protein MSAs |
| `--run_graft_phenotype` | `false` | ASR/selection/ranking; runs before the AF handoff |
| `--graft_matrix` | `data/graft_outcomes_11-09-2026.tsv` | Rootstock, scion, and TRUE/FALSE success |
| `--run_hyphy` | `false` | Selection tests (needs codon MSAs). `-profile conda` installs HyPhy from `conf/hyphy.yml` |
| `--run_asr` | `true` | IQ-TREE ancestral reconstruction |
| `--run_convergence` | `true` | PCOC and CSUBST; images are pulled if the programs are not on PATH |
| `--pcoc_image` | `carinerey/pcoc` | PCOC image used by the convergence step |
| `--csubst_image` | `quay.io/biocontainers/csubst:1.14.4--py311h04f3980_0` | CSUBST image used by the convergence step |
| `--prepare_alphafold` | `true` | Write the AF handoff as the last CPU step |
| `--af_min_priority_score` | `0` | Min phenotype score to stage an OG |
| `--af_max_ogs` | `0` | Cap OGs sent to GPU (`0` = no cap) |
| `--orthofinder_threads` | `16` | OrthoFinder CPUs |
| `--eggnog_cpus` | profile default | CPUs for `emapper.py` (`24` on `-profile hpc`) |
| `--eggnog_extra` | empty | Extra `emapper.py` args, e.g. faster Diamond search |

Database paths and any of the options above can be set in `conf/server.config` (copy `conf/server.config.example`). The matching flags are `--eggnog_data_dir`, `--interproscan_dir`, `--signalp_bin`, `--deeptmhmm_bin`, `--busco_download_path`, and `--usalign_bin`.

---

## Main outputs

| Path | Contents |
|------|----------|
| `01_raw_proteomes/` | Downloaded proteins, CDS, GFF, and `*.provenance.tsv` |
| `02_clean_proteomes/` | QC’d proteins + id maps |
| `02_clean_cds/` | QC’d CDS matching proteins |
| `02b_busco/` | BUSCO summaries + pass/fail gate logs |
| `03_orthofinder/` | Orthogroups, MSAs, species/gene trees |
| `05_selected_orthogroups/` | GOI orthogroups, member FASTA, `gene_id_crosswalk.tsv`, and `candidate_cds_check.tsv` when `--candidate_cds` is set |
| `05b_codon_alignments/` | Codon MSAs for selected OGs |
| `04_annotation/` | eggNOG / InterProScan / SignalP / DeepTMHMM |
| `06_percent_identity/` | Identity to the reference, plus an all-by-all matrix for each family |
| `07_tables/` | Polished family tables |
| `10_graft_phenotype/` | Traits, ASR, divergence tests, ranking (before the handoff) |
| `08_alphafold_inputs/` | GPU handoff package (written last; folder number kept) |
| `09_structure_comparison/` | RMSD, substitutions, superimposed PDBs (after GPU folding) |
| `11_curated/` | Family report (`index.html`), genome and BUSCO page (`genomes.html`), protein alignments (`protein_alignments.html`), codon alignments (`codon_alignments.html`), AlphaFold viewer (`structures.html`), and copies of the tables |

---

## Project layout

```text
GraftingGeneFamilies/
├── main.nf                          # Primary CPU pipeline
├── nextflow.config                  # Params + profiles
├── conf/
│   ├── species.yaml                 # Species, URLs, accessions, BUSCO lineage
│   ├── environment.yml              # Conda dependencies
│   ├── hyphy.yml                    # HyPhy env for the selection step
│   ├── server.config.example       # Copy to server.config for paths and options
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
│   ├── candidate_cds.example.fna    # Header example; nucleotides are placeholders
│   ├── grafting_candidates.txt      # Your GOI list on the server
│   ├── graft_outcomes.example.tsv
│   └── graft_outcomes_11-09-2026.tsv
├── docs/
│   ├── alphafold_plan.md
│   ├── graft_phenotype.md
│   └── structure_feature_mapping.md
└── assets/empty.tsv
```

---

## Server prerequisites

| Component | Notes |
|-----------|--------|
| Nextflow ≥ 23.04 | Installed on `curie` |
| Conda env | `conf/environment.yml` (OrthoFinder, BUSCO, IQ-TREE, Biopython, …). HyPhy is installed for the selection step from `conf/hyphy.yml` when `-profile conda` is on |
| eggNOG DB | `--eggnog_data_dir` |
| InterProScan | `--interproscan_dir` |
| SignalP 6 | Licensed binary; `--signalp_bin` |
| BUSCO lineages | Prefer shared `--busco_download_path` |
| DeepTMHMM | BioLib or local install |
| LocalColabFold | GPU server only |
| PCOC / CSUBST | Pulled during the convergence step (`carinerey/pcoc` and the CSUBST Biocontainer) when Docker or Singularity/Apptainer is available. The summary records a failed row if the pull cannot run |
| USalign / TMalign | Optional TM-score in structure comparison |

---

## Design choices (short)

- **Nextflow** for HPC resume, profiles, and process isolation  
- **ITAG / Niben community IDs** so GOI lists stay human-readable  
- **BUSCO gate** before OrthoFinder so bad annotations don’t pollute families  
- **Annotate selected OGs only** to save InterProScan/eggNOG time  
- **Optional characterized CDS check** so a gene-name match is confirmed against the sequence from the species where that gene was described
- **Codon MSAs from protein MSAs** so HyPhy uses the same homology hypothesis as OrthoFinder  
- **AlphaFold handoff last** so phenotype ranking can drop low-priority families before anything is staged for the GPU  
- **AlphaFold split CPU/GPU** so OrthoFinder and folding don’t fight for resources  
- **Phenotype tests control for phylogeny** (shift branches, partial Mantel, HyPhy foregrounds) rather than raw residue counting  

---

## Suggested next steps

1. Replace example GOIs / graft matrix with experimental data  
2. Point all DB/tool params at real paths on `curie`  
3. Complete a `-resume` run through tables + codon alignments  
4. Turn on `--run_graft_phenotype true` once the matrix is real  
5. After the CPU pipeline finishes, sync the AlphaFold handoff (its last output) to the GPU server and fold prioritized pairs  
6. Iterate on candidates that are (a) on phenotype-shift branches, (b) significant beyond species distance, and (c) structurally plausible in AF comparisons  

Questions or design changes: adjust `conf/species.yaml` and `nextflow.config` first — most behavior is parameterized there.
