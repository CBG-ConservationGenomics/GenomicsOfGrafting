process RUN_ORTHOFINDER {
    label 'process_highmem'
    publishDir "${params.outdir}/03_orthofinder", mode: 'copy'

    input:
    path proteomes

    output:
    path "OrthoFinder/Results_*", emit: results
    path "OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv", emit: orthogroups
    path "OrthoFinder/Results_*/Orthogroups/Orthogroups.txt", optional: true, emit: orthogroups_txt
    path "OrthoFinder/Results_*/MultipleSequenceAlignments", emit: alignments
    path "OrthoFinder/Results_*/Species_Tree/SpeciesTree_rooted.txt", emit: species_tree
    path "OrthoFinder/Results_*/Gene_Trees", optional: true, emit: gene_trees
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    mkdir -p proteome_input
    for f in ${proteomes}; do
        cp "\$f" proteome_input/
    done

    orthofinder \\
        -f proteome_input \\
        -t ${task.cpus} \\
        -a ${task.cpus} \\
        -M msa \\
        ${params.orthofinder_extra} \\
        -o OrthoFinder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder 2>&1 | head -n1 | sed 's/[^0-9.]*//')
    END_VERSIONS
    """
}
