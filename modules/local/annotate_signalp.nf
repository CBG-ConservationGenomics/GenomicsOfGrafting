process ANNOTATE_SIGNALP {
    label 'process_medium'
    publishDir "${params.outdir}/04_annotation/signalp", mode: 'copy'

    input:
    path fasta

    output:
    path "signalp.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    params.run_signalp

    script:
    """
    set -euo pipefail

    # SignalP 6.x typical CLI; adjust on your server if needed
    ${params.signalp_bin} \\
        --fastafile ${fasta} \\
        --output_dir signalp_out \\
        --format txt \\
        --organism euk \\
        --mode fast \\
        || ${params.signalp_bin} -fasta ${fasta} -stdout > signalp_out/prediction_results.txt

    # Normalize to a simple TSV if possible
    if [ -f signalp_out/prediction_results.txt ]; then
        cp signalp_out/prediction_results.txt signalp.tsv
    elif [ -f signalp_out/output.gff3 ]; then
        cp signalp_out/output.gff3 signalp.tsv
    else
        find signalp_out -type f | head
        # fallback: concatenate any tabular outputs
        cat signalp_out/* > signalp.tsv || touch signalp.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(${params.signalp_bin} --version 2>&1 | head -n1 || echo 'NA')
    END_VERSIONS
    """
}
