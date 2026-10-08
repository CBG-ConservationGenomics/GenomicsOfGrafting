process BUILD_CODON_ALIGNMENTS {
    label 'process_medium'
    publishDir "${params.outdir}/05b_codon_alignments", mode: 'copy'

    input:
    path og_list
    path alignments_dir
    path cds_fastas

    output:
    path "codon_alignments", emit: dir
    path "codon_alignment_report.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p codon_alignments

    build_codon_alignments.py \\
        --og-list ${og_list} \\
        --alignments-dir ${alignments_dir} \\
        --cds-fastas ${cds_fastas} \\
        --outdir codon_alignments \\
        --out-report codon_alignment_report.tsv \\
        --min-sequences ${params.codon_min_sequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
