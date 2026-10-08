process PREPARE_ALPHAFOLD_INPUTS {
    label 'process_low'
    publishDir "${params.outdir}/08_alphafold_inputs", mode: 'copy'

    input:
    path selected_ogs
    path identity_summary
    path identity_pairwise
    path proteins_fasta

    output:
    path "alphafold_inputs", emit: handoff
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    prepare_alphafold_inputs.py \\
        --selected-ogs ${selected_ogs} \\
        --identity-summary ${identity_summary} \\
        --identity-pairwise ${identity_pairwise} \\
        --proteins-fasta ${proteins_fasta} \\
        --outdir alphafold_inputs \\
        --min-pct-id ${params.af_min_pct_id} \\
        --max-models-per-og ${params.af_max_models_per_og} \\
        --mode ${params.af_input_mode}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
