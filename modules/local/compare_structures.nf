process COMPARE_STRUCTURES {
    label 'process_medium'
    publishDir "${params.outdir}/09_structure_comparison", mode: 'copy'

    input:
    path pairs
    path models_dir
    path alignments_dir

    output:
    path "structure_comparison_summary.tsv", emit: summary
    path "substitutions.tsv", emit: substitutions
    path "superimposed", optional: true, emit: superimposed
    path "versions.yml", emit: versions

    script:
    def usalign = params.usalign_bin ? "--usalign-bin ${params.usalign_bin}" : "--usalign-bin ''"
    """
    set -euo pipefail

    compare_structures.py \\
        --pairs ${pairs} \\
        --models-dir ${models_dir} \\
        --alignments-dir ${alignments_dir} \\
        --outdir . \\
        --min-plddt ${params.af_min_plddt} \\
        ${params.af_write_superimposed ? '--write-superimposed' : '--no-write-superimposed'} \\
        ${usalign}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
