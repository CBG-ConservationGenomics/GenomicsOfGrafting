process SELECT_ORTHOGROUPS {
    label 'process_low'
    publishDir "${params.outdir}/05_selected_orthogroups", mode: 'copy'

    input:
    path orthogroups_tsv
    path genes_of_interest
    path id_maps

    output:
    path "selected_orthogroups.tsv", emit: tsv
    path "selected_orthogroups.txt", emit: list
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    # Concatenate per-species id maps (keep one header)
    head -n1 \$(ls ${id_maps} | head -n1) > all_id_maps.tsv
    for f in ${id_maps}; do
        tail -n +2 "\$f" >> all_id_maps.tsv
    done

    select_orthogroups.py \\
        --orthogroups ${orthogroups_tsv} \\
        --genes-of-interest ${genes_of_interest} \\
        --id-map all_id_maps.tsv \\
        --out-tsv selected_orthogroups.tsv \\
        --out-og-list selected_orthogroups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
