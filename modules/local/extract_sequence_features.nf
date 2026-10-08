process EXTRACT_SEQUENCE_FEATURES {
    label 'process_low'
    publishDir "${params.outdir}/09_structure_comparison/features", mode: 'copy'

    input:
    path proteins
    path signalp
    path interproscan
    path deeptmhmm
    path eggnog
    path rules
    path transferred_sites

    output:
    path "sequence_features.tsv", emit: features
    path "protein_classifications.tsv", emit: classifications
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    SP=""; IPS=""; TM=""; EGG=""; XFER=""
    [ -s "${signalp}" ] && SP="--signalp ${signalp}"
    [ -s "${interproscan}" ] && IPS="--interproscan ${interproscan}"
    [ -s "${deeptmhmm}" ] && TM="--deeptmhmm ${deeptmhmm}"
    [ -s "${eggnog}" ] && EGG="--eggnog ${eggnog}"
    [ -s "${transferred_sites}" ] && XFER="--transferred-sites ${transferred_sites}"

    extract_sequence_features.py \\
        --proteins ${proteins} \\
        --rules ${rules} \\
        \$SP \$IPS \$TM \$EGG \$XFER \\
        --out-tsv sequence_features.tsv \\
        --out-classifications protein_classifications.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
