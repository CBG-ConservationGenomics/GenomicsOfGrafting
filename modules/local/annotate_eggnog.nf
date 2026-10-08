process ANNOTATE_EGGNOG {
    label 'process_high'
    publishDir "${params.outdir}/04_annotation/eggnog", mode: 'copy'

    input:
    path fasta

    output:
    path "eggnog.emapper.annotations", emit: annotations
    path "versions.yml", emit: versions

    when:
    params.run_eggnog

    script:
    def db = params.eggnog_data_dir ? "--data_dir ${params.eggnog_data_dir}" : ""
    """
    set -euo pipefail

    emapper.py \\
        -i ${fasta} \\
        -o eggnog \\
        --cpu ${task.cpus} \\
        --itype proteins \\
        ${db} \\
        --override

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
