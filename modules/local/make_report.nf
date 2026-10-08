process MAKE_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/07_tables", mode: 'copy'

    input:
    path selected_ogs
    path identity_summary
    path identity_pairwise
    path eggnog
    path interproscan
    path signalp
    path deeptmhmm

    output:
    path "gene_families_of_interest.tsv", emit: table
    path "gene_families_pctid_by_species.tsv", emit: wide
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    # Optional annotation files may be placeholders
    EGGNOG_ARG=""
    IPS_ARG=""
    SIGNALP_ARG=""
    TM_ARG=""
    [ -s "${eggnog}" ] && EGGNOG_ARG="--eggnog ${eggnog}"
    [ -s "${interproscan}" ] && IPS_ARG="--interproscan ${interproscan}"
    [ -s "${signalp}" ] && SIGNALP_ARG="--signalp ${signalp}"
    [ -s "${deeptmhmm}" ] && TM_ARG="--deeptmhmm ${deeptmhmm}"

    make_family_table.py \\
        --selected-ogs ${selected_ogs} \\
        --identity-summary ${identity_summary} \\
        --identity-pairwise ${identity_pairwise} \\
        \$EGGNOG_ARG \$IPS_ARG \$SIGNALP_ARG \$TM_ARG \\
        --out-table gene_families_of_interest.tsv \\
        --out-wide-identity gene_families_pctid_by_species.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
