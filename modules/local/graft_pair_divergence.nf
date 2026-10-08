process GRAFT_PAIR_DIVERGENCE {
    label 'process_medium'
    publishDir "${params.outdir}/10_graft_phenotype/divergence", mode: 'copy'

    input:
    path pairs
    path og_list
    path alignments_dir

    output:
    path "graft_pair_divergence.tsv", emit: divergence
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    calc_graft_pair_divergence.py \\
        --pairs ${pairs} \\
        --og-list ${og_list} \\
        --alignments-dir ${alignments_dir} \\
        --out-tsv graft_pair_divergence.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
