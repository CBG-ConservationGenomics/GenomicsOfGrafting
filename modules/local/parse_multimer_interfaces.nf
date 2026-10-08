process PARSE_MULTIMER_INTERFACES {
    label 'process_low'
    publishDir "${params.outdir}/09_structure_comparison/multimer_interfaces", mode: 'copy'

    input:
    path models_dir
    path manifest

    output:
    path "multimer_interface_summary.tsv", emit: summary
    path "compatible_vs_incompatible_interface.tsv", emit: comparison
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    parse_multimer_interfaces.py \\
        --models-dir ${models_dir} \\
        --manifest ${manifest} \\
        --outdir . \\
        --distance-cutoff ${params.interface_distance_cutoff} \\
        --pae-cutoff ${params.interface_pae_cutoff} \\
        --min-plddt ${params.af_min_plddt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
