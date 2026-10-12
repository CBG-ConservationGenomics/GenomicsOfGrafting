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
    def extra = params.eggnog_extra?.trim() ? " ${params.eggnog_extra.trim()}" : ""
    """
    set -euo pipefail

    if [ ! -s "${fasta}" ]; then
        echo "ERROR: ${fasta} is empty, so eggNOG was not started. No orthogroup proteins were extracted." >&2
        exit 1
    fi

    emapper.py \\
        -i ${fasta} \\
        -o eggnog \\
        --cpu ${task.cpus} \\
        --itype proteins \\
        ${db}${extra} \\
        --override

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
