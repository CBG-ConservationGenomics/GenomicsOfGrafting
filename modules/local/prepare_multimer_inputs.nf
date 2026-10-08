process PREPARE_MULTIMER_INPUTS {
    label 'process_low'
    publishDir "${params.outdir}/08_alphafold_inputs", mode: 'copy'

    input:
    path pairs
    path classifications
    path proteins
    path selected_ogs

    output:
    path "multimer", emit: package
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    prepare_multimer_inputs.py \\
        --pairs ${pairs} \\
        --classifications ${classifications} \\
        --proteins ${proteins} \\
        --selected-ogs ${selected_ogs} \\
        --outdir multimer

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
