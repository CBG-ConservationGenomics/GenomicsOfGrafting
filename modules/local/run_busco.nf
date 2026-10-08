process RUN_BUSCO {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/02b_busco/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.busco_metrics.tsv"), emit: metrics
    tuple val(meta), path("${meta.id}.busco_status.txt"), emit: status
    path "busco_run", emit: busco_dir
    path "versions.yml", emit: versions

    script:
    def lineage = meta.busco_lineage ?: params.busco_lineage
    def download = params.busco_download_path ? "--download_path ${params.busco_download_path}" : ""
    def offline = params.busco_offline ? "--offline" : ""
    """
    set -euo pipefail

    busco \\
        -i ${fasta} \\
        -l ${lineage} \\
        -m proteins \\
        -o busco_run \\
        -c ${task.cpus} \\
        --force \\
        ${download} \\
        ${offline}

    SUMMARY=\$(find busco_run -name 'short_summary*.txt' | head -n1)
    if [ -z "\$SUMMARY" ]; then
        echo "ERROR: BUSCO short_summary not found for ${meta.id}" >&2
        find busco_run -type f | head -n50 >&2 || true
        exit 1
    fi

    parse_busco_summary.py \\
        --summary "\$SUMMARY" \\
        --species-id ${meta.id} \\
        --lineage ${lineage} \\
        --min-complete ${params.busco_min_complete} \\
        --out-tsv ${meta.id}.busco_metrics.tsv \\
        --out-status ${meta.id}.busco_status.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
