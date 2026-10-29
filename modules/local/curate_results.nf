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
    path protein_alignments
    path candidate_check
    path species_tree

    output:
    path "index.html", emit: html
    path "genomes.html", emit: genomes
    path "protein_alignments.html", emit: proteins
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

    CHECK_ARG=""
    if [ -s "${candidate_check}" ] && grep -q query_gene "${candidate_check}" 2>/dev/null; then
        CHECK_ARG="--candidate-check ${candidate_check}"
    fi

    TREE_ARG=""
    if [ -s "${species_tree}" ]; then
        TREE_ARG="--species-tree ${species_tree}"
    fi

    curate_results.py \\
        --families ${families} \\
        --wide-identity ${wide} \\
        \$RANK_ARG \\
        \$AF_ARG \\
        \$CHECK_ARG \\
        --provenance ${provenance_files} \\
        --identity-matrix-pairs ${matrix_pairs} \\
        --identity-matrix-labels ${matrix_labels} \\
        \$TREE_ARG \\
        --runtime-seconds ${runtimeSeconds} \\
        --outdir .

    write_summary_pages.py \\
        --provenance ${provenance_files} \\
        --busco ${busco_summary} \\
        --codon-report ${codon_report} \\
        --codon-dir ${codon_dir} \\
        --families ${families} \\
        --protein-alignments ${protein_alignments} \\
        --outdir .
    # Family table cells are the closest protein to the reference, not the species average.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
