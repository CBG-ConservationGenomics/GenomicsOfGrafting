#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * GPU-server AlphaFold / ColabFold workflow
 *
 * Run this ON the GPU machine after copying results/08_alphafold_inputs/alphafold_inputs:
 *
 *   nextflow run workflows/alphafold.nf -profile gpu \
 *     --alphafold_inputs /data/handoff/alphafold_inputs \
 *     --outdir alphafold_results
 */

params.alphafold_inputs = null
params.outdir           = 'alphafold_results'
params.af_tool          = 'colabfold'   // colabfold for now
params.af_num_models    = 1
params.af_num_recycle   = 3
params.af_use_amber     = false

if (!params.alphafold_inputs) {
    error "Set --alphafold_inputs to the handoff directory (…/08_alphafold_inputs/alphafold_inputs)"
}

process RUN_COLABFOLD {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/models", mode: 'copy'
    label 'gpu'

    input:
    path fasta

    output:
    path "${fasta.baseName}_out", emit: model_dir

    script:
    def amber = params.af_use_amber ? '--amber' : ''
    """
    set -euo pipefail
    mkdir -p ${fasta.baseName}_out

    if command -v colabfold_batch >/dev/null 2>&1; then
        colabfold_batch \\
            ${fasta} \\
            ${fasta.baseName}_out \\
            --num-models ${params.af_num_models} \\
            --num-recycle ${params.af_num_recycle} \\
            ${amber}
    else
        echo "ERROR: colabfold_batch not found on PATH." >&2
        echo "Install LocalColabFold on this GPU server and retry." >&2
        exit 1
    fi
    """
}

process SUMMARIZE_AF {
    publishDir "${params.outdir}", mode: 'copy'
    label 'process_low'

    input:
    path model_dirs
    path manifest

    output:
    path "alphafold_summary.tsv"

    script:
    """
    set -euo pipefail
    summarize_alphafold.py \\
        --manifest ${manifest} \\
        --model-dirs ${model_dirs} \\
        --out-tsv alphafold_summary.tsv
    """
}

workflow {
    ch_fastas   = Channel.fromPath("${params.alphafold_inputs}/sequences/*.fasta", checkIfExists: true)
    ch_manifest = Channel.fromPath("${params.alphafold_inputs}/manifest.tsv", checkIfExists: true)

    if (params.af_tool != 'colabfold') {
        error "params.af_tool=${params.af_tool} not implemented yet; use --af_tool colabfold"
    }

    RUN_COLABFOLD(ch_fastas)
    SUMMARIZE_AF(RUN_COLABFOLD.out.model_dir.collect(), ch_manifest)
}

workflow.onComplete {
    log.info "AlphaFold workflow complete → ${params.outdir}"
}
