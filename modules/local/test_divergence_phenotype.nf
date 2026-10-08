process TEST_DIVERGENCE_PHENOTYPE {
    label 'process_low'
    publishDir "${params.outdir}/10_graft_phenotype/divergence", mode: 'copy'

    input:
    path pair_divergence
    path species_distance

    output:
    path "divergence_vs_phenotype.tsv", emit: tests
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    test_divergence_vs_phenotype.py \\
        --pair-divergence ${pair_divergence} \\
        --species-distance ${species_distance} \\
        --out-tsv divergence_vs_phenotype.tsv \\
        --n-perm ${params.graft_mantel_perms}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
