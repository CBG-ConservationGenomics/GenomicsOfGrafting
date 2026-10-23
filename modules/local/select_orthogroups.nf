process SELECT_ORTHOGROUPS {
    label 'process_low'
    publishDir "${params.outdir}/05_selected_orthogroups", mode: 'copy'

    input:
    path orthofinder_results
    path genes_of_interest
    path id_maps
    path gffs

    output:
    path "selected_orthogroups.tsv", emit: tsv
    path "selected_orthogroups.txt", emit: list
    path "gene_id_crosswalk.tsv", emit: crosswalk
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    # Concatenate per-species id maps (keep one header)
    head -n1 \$(ls ${id_maps} | head -n1) > all_id_maps.tsv
    for f in ${id_maps}; do
        tail -n +2 "\$f" >> all_id_maps.tsv
    done

    # BAC gene names stay whole, and Note cross-references are not locus ids.
    select_orthogroups.py \\
        --orthogroups ${orthofinder_results} \\
        --genes-of-interest ${genes_of_interest} \\
        --id-map all_id_maps.tsv \\
        --gff ${gffs} \\
        --out-tsv selected_orthogroups.tsv \\
        --out-og-list selected_orthogroups.txt \\
        --out-crosswalk gene_id_crosswalk.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
