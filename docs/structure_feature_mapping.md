# Mapping functional domains onto AlphaFold structures

Connect sequence annotations → residue coordinates → AlphaFold models for grafting-relevant protein classes.

## Core idea

AlphaFold models are numbered by the **input sequence** (residue 1…N).  
Annotations on that same sequence map 1:1 onto the structure.  
Across species, transfer sites via **OrthoFinder MSA columns**, never by copying raw residue numbers.

```text
InterProScan / SignalP / DeepTMHMM / curated sites / selection
        │
        ▼
 sequence_features.tsv  (+ transferred_sites.tsv via MSA)
        │
        ▼
 AlphaFold PDB/CIF       (B-factor = pLDDT)
        │
        ├── monomers: enzymes / defense surface / transporters
        └── multimers: interface contacts + PAE
```

Run the integrated workflow:

```bash
nextflow run workflows/map_structure_features.nf -profile hpc,conda \
  --proteins results/05_selected_orthogroups/selected_og_proteins.faa \
  --signalp results/04_annotation/signalp/signalp.tsv \
  --interproscan results/04_annotation/interproscan/interproscan.tsv \
  --deeptmhmm results/04_annotation/deeptmhmm/deeptmhmm.tsv \
  --eggnog results/04_annotation/eggnog/eggnog.emapper.annotations \
  --models_dir alphafold_results/models \
  --substitutions results/09_structure_comparison/substitutions.tsv \
  --alignments_dir results/03_orthofinder/.../MultipleSequenceAlignments \
  --selected_ogs results/05_selected_orthogroups/selected_orthogroups.tsv \
  --multimer_models_dir alphafold_multimer/models \
  --multimer_manifest results/08_alphafold_inputs/multimer/manifest.tsv \
  --outdir results
```

---

## 1. Multimer interface PAE parser

**Goal:** For each scion-ligand × rootstock-receptor ColabFold multimer job, report contact count, mean/min interface PAE, interface pLDDT, iptm/ptm; then compare compatible vs incompatible grafts.

### Inputs

| File | Role |
|------|------|
| `multimer/manifest.tsv` | From `prepare_multimer_inputs.py` (`job_id`, `ligand_id`, `receptor_id`, `success`, …) |
| ColabFold output dirs | Ranked `.pdb` + `*_scores_rank_*.json` containing `pae` / `plddt` / `iptm` |

### Standalone

```bash
parse_multimer_interfaces.py \
  --models-dir alphafold_multimer/models \
  --manifest results/08_alphafold_inputs/multimer/manifest.tsv \
  --outdir results/09_structure_comparison/multimer_interfaces \
  --distance-cutoff 5.0 \
  --pae-cutoff 15.0 \
  --min-plddt 50
```

### What it does

1. Finds best-ranked PDB per `job_id` and matching scores JSON.  
2. Enumerates heavy-atom contacts between chain A and B (≤ `--distance-cutoff` Å).  
3. Looks up PAE at absolute residue indices (concatenated chain order, matching AF/ColabFold).  
4. Counts **confident** contacts: PAE ≤ cutoff and mean pLDDT ≥ min.  
5. Writes:
   - `multimer_interface_summary.tsv` — per job  
   - `compatible_vs_incompatible_interface.tsv` — mean metrics by `success=yes` vs `no`

### Interpretation

Lower mean interface PAE and more confident contacts in compatible pairs is a **testable structural hypothesis**, not proof of graft compatibility.

Upstream job prep: `prepare_multimer_inputs.py` + GPU `workflows/alphafold.nf` multimer mode.

---

## 2. Surface accessibility for defense proteins

**Goal:** Flag solvent-exposed residues on AF models; highlight scion–rootstock substitutions on the surface (combine with HyPhy diversifying sites).

### Method order

1. **Biopython `ShrakeRupley`** (preferred; biopython ≥1.80 in `environment.yml`)  
2. **FreeSASA** Python API if installed (`pip install freesasa` / conda)  
3. **CA/CB neighbor density** fallback (exposure proxy, not true RSA)

RSA = residue SASA / max ASA (Tien 2013). Default surface: RSA ≥ `0.25`.

### Standalone

```bash
# After extract_sequence_features → protein_classifications.tsv
compute_surface_rsa.py \
  --models-dir alphafold_results/models \
  --classifications results/09_structure_comparison/features/protein_classifications.tsv \
  --substitutions results/09_structure_comparison/substitutions.tsv \
  --outdir results/09_structure_comparison/surface \
  --rsa-surface-cutoff 0.25
```

If classifications include a `defense` class (from `conf/functional_feature_rules.yaml` keywords), only those proteins are processed. Otherwise all IDs in the classifications file are used.

### Outputs

| File | Contents |
|------|----------|
| `residue_surface.tsv` | protein, resnum, aa, score, is_surface, method |
| `substitutions_on_surface.tsv` | substitutions with `ref_is_surface` / `orth_is_surface` |

**Priority for defense:** surface ∩ selected (MEME/FEL) ∩ graft-partner difference ≫ buried conserved core.

---

## 3. Curated active-site / pore residues via OrthoFinder MSA

**Goal:** Curate catalytic or pore residues on **one reference sequence**, then map every ortholog to its own residue numbers through the OG MSA.

### Why not raw numbers across species?

Tomato residue 198 is not Arabidopsis residue 198. Gaps and indels shift numbering. The MSA column is the homology link.

### Step A — curate in YAML

Edit `conf/functional_feature_rules.yaml`:

```yaml
curated_sites:
  - orthogroup: OG0001234          # OrthoFinder OG id
    reference_id: "Sly|Solanum_lycopersicum|Solyc03g116610.3.1"
    feature: active_site           # or pore_site
    residues: [198, 201, 310]      # 1-based on the REFERENCE ungapped protein only
    note: "GH9 catalytic residues from literature on this tomato ID"
```

### Step B — transfer

```bash
transfer_curated_sites.py \
  --rules conf/functional_feature_rules.yaml \
  --alignments-dir results/03_orthofinder/.../MultipleSequenceAlignments \
  --selected-ogs results/05_selected_orthogroups/selected_orthogroups.tsv \
  --out-tsv results/09_structure_comparison/features/transferred_sites.tsv
```

Logic per curated residue:

```text
ref ungapped position → MSA column → each ortholog’s ungapped position
(status = gap_in_target if that column is a gap)
```

### Step C — merge into structure mapping

`extract_sequence_features.py --transferred-sites transferred_sites.tsv` appends point features (`source=curated_msa_transfer`).  
`map_features_to_structure.py` then flags substitutions inside those features on AF models.

The workflow runs A→C automatically when `--alignments_dir` is set.

### Status codes in `transferred_sites.tsv`

| status | Meaning |
|--------|---------|
| `ok` | Mapped to a non-gap residue in the target |
| `gap_in_target` | Homologous column is a gap in that ortholog |
| `reference_not_in_msa` | Fix `reference_id` to match OrthoFinder headers |
| `alignment_missing` | OG MSA not found under `--alignments-dir` |
| `ref_residue_out_of_range` | Curated number longer than reference sequence |

---

## Workflow by protein class (summary)

### A. Cell-wall enzymes

SignalP + InterPro + **MSA-transferred active sites** + glyco sequons → flag substitutions in catalytic/binding features.

### B. Ligand–receptor

`prepare_multimer_inputs.py` → ColabFold multimer → **`parse_multimer_interfaces.py`**.

### C. Defense

Codon MSA + HyPhy + **`compute_surface_rsa.py`** + substitutions on surface.

### D. Transporters / PD

DeepTMHMM + **MSA-transferred pore sites** + AF monomer.

---

## Scripts & modules

| Script / module | Role |
|-----------------|------|
| `conf/functional_feature_rules.yaml` | Class keywords + curated site lists |
| `bin/transfer_curated_sites.py` | Reference residues → all OG members via MSA |
| `bin/extract_sequence_features.py` | Unify annotations (+ optional transferred sites) |
| `bin/map_features_to_structure.py` | Features × PDB × substitutions |
| `bin/compute_surface_rsa.py` | RSA / surface for defense |
| `bin/parse_multimer_interfaces.py` | Interface contacts + PAE |
| `bin/prepare_multimer_inputs.py` | Cross-partner multimer FASTAs |
| `workflows/map_structure_features.nf` | Orchestrates the above |

Outputs under `results/09_structure_comparison/{features,surface,multimer_interfaces}/`.

---

## Practical rules

1. AF input sequence = annotated renamed FASTA (`Sly|…|id`).  
2. Cross-species site homology = MSA column, not shared residue index.  
3. Down-weight pLDDT &lt; 50 for mechanistic claims.  
4. Signal peptides are often cleaved/disordered — prioritize mature domains unless trafficking is the hypothesis.
