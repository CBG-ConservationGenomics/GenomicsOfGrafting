#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * CPU-side structure comparison after AlphaFold/ColabFold.
 *
 *   nextflow run workflows/compare_structures.nf -profile hpc,conda \
 *     --alphafold_inputs results/08_alphafold_inputs/alphafold_inputs \
 *     --alphafold_models alphafold_results/models \
 *     --alignments_dir results/03_orthofinder/.../MultipleSequenceAlignments \
 *     --outdir results
 */

include { COMPARE_STRUCTURES } from '../modules/local/compare_structures'

params.alphafold_inputs = null
params.alphafold_models = null
params.alignments_dir   = null
params.outdir           = 'results'
params.af_min_plddt     = 50.0
params.af_write_superimposed = true
params.usalign_bin      = 'USalign'   // set to '' to skip TM-score

if (!params.alphafold_inputs) {
    error "Set --alphafold_inputs (directory containing pairs.tsv)"
}
if (!params.alphafold_models) {
    error "Set --alphafold_models (ColabFold/AF output directory)"
}
if (!params.alignments_dir) {
    error "Set --alignments_dir (OrthoFinder MultipleSequenceAlignments)"
}

workflow {
    ch_pairs  = Channel.fromPath("${params.alphafold_inputs}/pairs.tsv", checkIfExists: true)
    ch_models = Channel.fromPath(params.alphafold_models, checkIfExists: true)
    ch_aln    = Channel.fromPath(params.alignments_dir, checkIfExists: true)

    COMPARE_STRUCTURES(ch_pairs, ch_models, ch_aln)
}

workflow.onComplete {
    log.info """
    Structure comparison complete.
      Summary : ${params.outdir}/09_structure_comparison/structure_comparison_summary.tsv
      Subs    : ${params.outdir}/09_structure_comparison/substitutions.tsv
    """.stripIndent()
}
