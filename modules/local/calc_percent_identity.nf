process CALC_PERCENT_IDENTITY {
    label 'process_medium'
    publishDir "${params.outdir}/06_percent_identity", mode: 'copy'

    input:
    path og_list
    path selected_tsv
    path alignments_dir

    output:
    path "pairwise_percent_identity.tsv", emit: pairwise
    path "summary_percent_identity.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    calc_percent_identity.py \\
        --og-list ${og_list} \\
        --alignments-dir ${alignments_dir} \\
        --selected-ogs-tsv ${selected_tsv} \\
        --reference-species-id ${params.reference_species_id} \\
        --out-pairwise pairwise_percent_identity.tsv \\
        --out-summary summary_percent_identity.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
