process COLLECT_BUSCO_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/02b_busco", mode: 'copy'

    input:
    path metrics_files

    output:
    path "busco_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail

    head -n1 \$(ls ${metrics_files} | head -n1) > busco_summary.tsv
    for f in ${metrics_files}; do
        tail -n +2 "\$f" >> busco_summary.tsv
    done
    """
}
