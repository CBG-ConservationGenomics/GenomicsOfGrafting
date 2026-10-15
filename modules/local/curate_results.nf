process CURATE_RESULTS {
    label 'process_low'
    publishDir "${params.outdir}/11_curated", mode: 'copy'

    input:
    path families
    path wide
    path ranking
    path handoff
    path provenance_files
    path matrix_pairs
    path matrix_labels
    path busco_summary
    path codon_report
    path codon_dir

    output:
    path "index.html", emit: html
    path "genomes.html", emit: genomes
    path "codon_alignments.html", emit: codons
    path "tables", emit: tables
    path "README.txt", emit: readme
    path "versions.yml", emit: versions

    script:
    def runtimeSeconds = java.time.Duration.between(workflow.start.toInstant(), java.time.Instant.now()).toSeconds()
    """
    set -euo pipefail

    RANK_ARG=""
    AF_ARG=""
    if [ -s "${ranking}" ] && grep -q Orthogroup "${ranking}" 2>/dev/null; then
        RANK_ARG="--ranking ${ranking}"
    fi
    if [ -d "${handoff}" ] && [ -f "${handoff}/og_priority.tsv" ]; then
        AF_ARG="--alphafold-dir ${handoff}"
    fi

    curate_results.py \\
        --families ${families} \\
        --wide-identity ${wide} \\
        \$RANK_ARG \\
        \$AF_ARG \\
        --provenance ${provenance_files} \\
        --identity-matrix-pairs ${matrix_pairs} \\
        --identity-matrix-labels ${matrix_labels} \\
        --runtime-seconds ${runtimeSeconds} \\
        --outdir .

    write_summary_pages.py \\
        --provenance ${provenance_files} \\
        --busco ${busco_summary} \\
        --codon-report ${codon_report} \\
        --codon-dir ${codon_dir} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
