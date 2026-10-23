process CHECK_CANDIDATE_CDS {
    label 'process_low'
    publishDir "${params.outdir}/05_selected_orthogroups", mode: 'copy', pattern: "candidate_cds_check.*"

    input:
    path crosswalk
    path proteins
    path cds_files
    path reference_cds

    output:
    path "candidate_cds_check.tsv", emit: report
    path "candidate_cds_check.summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    def require_flag = params.candidate_cds_require_all ? "--require-all" : ""
    """
    set -euo pipefail

    check_candidate_cds.py \\
        --crosswalk ${crosswalk} \\
        --proteins ${proteins} \\
        --cds ${cds_files} \\
        --reference-cds ${reference_cds} \\
        --min-protein-identity ${params.candidate_cds_min_protein_identity} \\
        --min-cds-identity ${params.candidate_cds_min_cds_identity} \\
        ${require_flag} \\
        --out candidate_cds_check.tsv \\
        --summary candidate_cds_check.summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
