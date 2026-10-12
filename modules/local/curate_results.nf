process CURATE_RESULTS {
    label 'process_low'
    publishDir "${params.outdir}/11_curated", mode: 'copy'

    input:
    path families
    path wide
    path ranking
    path handoff

    output:
    path "index.html", emit: html
    path "tables", emit: tables
    path "README.txt", emit: readme
    path "versions.yml", emit: versions

    script:
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
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
