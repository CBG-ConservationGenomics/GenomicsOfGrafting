process MAP_FEATURES_TO_STRUCTURE {
    label 'process_low'
    publishDir "${params.outdir}/09_structure_comparison/features", mode: 'copy'

    input:
    path features
    path models_dir
    path substitutions

    output:
    path "structure_features.tsv", emit: structure_features
    path "substitutions_in_features.tsv", emit: subs_in_features
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    SUB=""
    [ -s "${substitutions}" ] && SUB="--substitutions ${substitutions}"

    map_features_to_structure.py \\
        --features ${features} \\
        --models-dir ${models_dir} \\
        \$SUB \\
        --outdir . \\
        --min-plddt ${params.af_min_plddt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
