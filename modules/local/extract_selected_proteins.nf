process EXTRACT_SELECTED_PROTEINS {
    label 'process_low'
    publishDir "${params.outdir}/05_selected_orthogroups", mode: 'copy'

    input:
    path selected_tsv
    path proteomes

    output:
    path "selected_og_proteins.faa", emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    extract_selected_proteins.py \\
        --selected-ogs ${selected_tsv} \\
        --proteomes ${proteomes} \\
        --out-fasta selected_og_proteins.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
