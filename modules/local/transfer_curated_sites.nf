process TRANSFER_CURATED_SITES {
    label 'process_low'
    publishDir "${params.outdir}/09_structure_comparison/features", mode: 'copy'

    input:
    path rules
    path alignments_dir
    path selected_ogs

    output:
    path "transferred_sites.tsv", emit: transferred
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    OG=""
    [ -s "${selected_ogs}" ] && OG="--selected-ogs ${selected_ogs}"

    transfer_curated_sites.py \\
        --rules ${rules} \\
        --alignments-dir ${alignments_dir} \\
        \$OG \\
        --out-tsv transferred_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
