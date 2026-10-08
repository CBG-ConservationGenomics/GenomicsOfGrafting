process ANNOTATE_INTERPROSCAN {
    label 'process_high'
    publishDir "${params.outdir}/04_annotation/interproscan", mode: 'copy'

    input:
    path fasta

    output:
    path "interproscan.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    params.run_interproscan

    script:
    def ips = params.interproscan_dir ? "${params.interproscan_dir}/interproscan.sh" : "interproscan.sh"
    """
    set -euo pipefail

    # InterProScan dislikes '*' stop codons in some builds — strip them
    sed 's/\\*//g' ${fasta} > proteins_nostop.faa

    ${ips} \\
        -i proteins_nostop.faa \\
        -f tsv \\
        -o interproscan.tsv \\
        -cpu ${task.cpus} \\
        -goterms \\
        -iprlookup \\
        -pa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(${ips} --version 2>&1 | head -n1 || echo 'NA')
    END_VERSIONS
    """
}
