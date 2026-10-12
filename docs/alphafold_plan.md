# AlphaFold integration (two-stage + CPU comparison)

## Recommended design

```text
CPU server (main.nf)
  … → gene family tables
  → graft phenotype (ASR / divergence / optional HyPhy)
  → 08_alphafold_inputs/   handoff ranked by phenotype priority
        ↓ rsync
GPU server (workflows/alphafold.nf)
  → ColabFold models (fold high fold_rank first)
        ↓ rsync models back
CPU server (workflows/compare_structures.nf)
  → RMSD, substitutions mapped to structure, optional TM-score
```

## 1. CPU handoff (`main.nf`)

Runs **after** `--run_graft_phenotype` when that flag is on, so
`graft_candidate_ranking.tsv` can prioritize which OGs to stage.

Writes `results/08_alphafold_inputs/alphafold_inputs/`:

| File | Purpose |
|------|---------|
| `sequences/*.fasta` | One protein per model |
| `manifest.tsv` | What to fold (sorted by `fold_rank`) |
| `og_priority.tsv` | OG-level phenotype scores used for ranking |
| `pairs.tsv` | Reference ↔ ortholog pairs for comparison |
| `README.txt` | GPU instructions |

Params: `--prepare_alphafold`, `--af_input_mode`, `--af_min_pct_id`,
`--af_max_models_per_og`, `--af_min_priority_score`, `--af_max_ogs`.

Priority score (from `rank_graft_candidates.py`) rewards ASR substitutions on
phenotype-shift branches, divergence beyond relatedness, HyPhy hits, and
convergence. Use `--af_max_ogs 20` or `--af_min_priority_score 5` to shrink the
GPU queue.

## 2. GPU folding (`workflows/alphafold.nf`)

```bash
rsync -av results/08_alphafold_inputs/alphafold_inputs/ gpu:/data/handoff/alphafold_inputs/
nextflow run workflows/alphafold.nf -profile gpu \
  --alphafold_inputs /data/handoff/alphafold_inputs \
  --outdir alphafold_results
```

Prefer **LocalColabFold**. Tune `conf/gpu.config` for your GPU queue/container.

## 3. CPU structure comparison (`workflows/compare_structures.nf`)

After models are back on the CPU server:

```bash
nextflow run workflows/compare_structures.nf -profile hpc,conda \
  --alphafold_inputs results/08_alphafold_inputs/alphafold_inputs \
  --alphafold_models alphafold_results/models \
  --alignments_dir results/03_orthofinder/OrthoFinder/Results_*/MultipleSequenceAlignments \
  --outdir results
```

Or call the script directly:

```bash
bin/compare_structures.py \
  --pairs results/08_alphafold_inputs/alphafold_inputs/pairs.tsv \
  --models-dir alphafold_results/models \
  --alignments-dir results/03_orthofinder/.../MultipleSequenceAlignments \
  --outdir results/09_structure_comparison \
  --min-plddt 50
```

### What it does

1. Reads each reference↔ortholog row from `pairs.tsv`
2. Finds the best-ranked PDB/CIF for each partner under the models directory
3. Uses the OrthoFinder MSA to define residue correspondence
4. Lists amino-acid substitutions with structure residue numbers + pLDDT (from B-factors)
5. Superposes CA atoms (Biopython) → `ca_rmsd` and optional `superimposed/*.pdb`
6. If `USalign` or `TMalign` is on PATH, records TM-scores

### Outputs (`results/09_structure_comparison/`)

| File | Contents |
|------|----------|
| `structure_comparison_summary.tsv` | Per-pair RMSD, n substitutions, TM-score |
| `substitutions.tsv` | Each AA change with MSA column, resnums, pLDDT, `low_confidence` flag |
| `superimposed/` | Ortholog PDBs fitted onto the reference |

### Params

- `--af_min_plddt` (default 50) — flag substitutions where either side is below this
- `--af_write_superimposed true|false`
- `--usalign_bin USalign` — set `''` to skip TM-score

## Functional feature mapping & multimers

See [`structure_feature_mapping.md`](structure_feature_mapping.md).

- `workflows/map_structure_features.nf` — SignalP/InterPro/TM/glyco → structure residues; MSA-transferred curated sites; defense surface RSA; multimer interface PAE  
- `bin/prepare_multimer_inputs.py` — scion ligand × rootstock receptor FASTAs for ColabFold multimer  
- `bin/parse_multimer_interfaces.py` — contacts + mean/min interface PAE; compatible vs incompatible summary  
- `bin/compute_surface_rsa.py` — Biopython/FreeSASA RSA (defense surface)  
- `bin/transfer_curated_sites.py` — curated active-site/pore residues → orthologs via OrthoFinder MSA  


## Later enhancements

- Pocket/interface proximity flags from multimer contacts  
- ChimeraX / PyMOL session scripts coloring feature-overlapping substitutions
