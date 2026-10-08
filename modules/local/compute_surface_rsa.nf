process COMPUTE_SURFACE_RSA {
    label 'process_low'
    publishDir "${params.outdir}/09_structure_comparison/surface", mode: 'copy'

    input:
    path models_dir
    path classifications
    path substitutions

    output:
    path "residue_surface.tsv", emit: surface
    path "substitutions_on_surface.tsv", emit: surface_subs
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    SUB=""
    [ -s "${substitutions}" ] && SUB="--substitutions ${substitutions}"

    compute_surface_rsa.py \\
        --models-dir ${models_dir} \\
        --classifications ${classifications} \\
        \$SUB \\
        --outdir . \\
        --rsa-surface-cutoff ${params.rsa_surface_cutoff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
