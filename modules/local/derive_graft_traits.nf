process DERIVE_GRAFT_TRAITS {
    label 'process_low'
    publishDir "${params.outdir}/10_graft_phenotype", mode: 'copy'

    input:
    path matrix
    path species_tree

    output:
    path "traits/pair_outcomes.tsv", emit: pairs
    path "traits/lineage_traits.tsv", emit: lineage
    path "traits/phenotype_shift_branches.tsv", emit: shift_branches
    path "traits/species_tree.annotated.nwk", emit: tree
    path "traits/species_distance.tsv", emit: sp_distance
    path "traits/foreground_tips.txt", emit: foreground
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p traits

    derive_graft_traits.py \\
        --matrix ${matrix} \\
        --species-tree ${species_tree} \\
        --outdir traits \\
        --broad-threshold ${params.graft_broad_threshold} \\
        --foreground-trait ${params.graft_foreground_trait}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
