process RUN_CONVERGENCE {
    tag "${og_id}"
    label 'process_medium'
    publishDir "${params.outdir}/10_graft_phenotype/convergence/${og_id}", mode: 'copy'

    input:
    tuple val(og_id), path(alignment), path(tree)
    path foreground_tips

    output:
    tuple val(og_id), path("${og_id}.convergence_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    params.run_convergence

    script:
    """
    set -euo pipefail

    echo -e "Orthogroup\\tmethod\\tsites_or_stat\\tstatus\\tnotes" > ${og_id}.convergence_summary.tsv

    prep_og_alignment.py \\
        --alignment ${alignment} \\
        --out-fasta ${og_id}.species.faa \\
        --out-species-map ${og_id}.species_map.tsv

    # --- PCOC (optional; requires pcoc_detect on PATH) ---
    if command -v pcoc_detect >/dev/null 2>&1 || command -v pcoc >/dev/null 2>&1; then
        # Write convergent tip set for PCOC
        cp ${foreground_tips} ${og_id}.convergent_leaves.txt
        if pcoc_detect -m ${og_id}.species.faa -t ${tree} -l ${og_id}.convergent_leaves.txt -o pcoc_out \\
            > pcoc.log 2>&1; then
            echo -e "${og_id}\\tPCOC\\tsee_pcoc_out\\tok\\t" >> ${og_id}.convergence_summary.tsv
        else
            echo -e "${og_id}\\tPCOC\\t\\tfailed\\tsee_pcoc.log" >> ${og_id}.convergence_summary.tsv
        fi
    else
        echo -e "${og_id}\\tPCOC\\t\\tskipped\\tpcoc_detect not installed" >> ${og_id}.convergence_summary.tsv
    fi

    # --- CSUBST (optional; requires csubst on PATH) ---
    if command -v csubst >/dev/null 2>&1; then
        if csubst analyze --alignment ${og_id}.species.faa --tree ${tree} \\
            --foreground ${foreground_tips} --outdir csubst_out > csubst.log 2>&1; then
            echo -e "${og_id}\\tCSUBST\\tsee_csubst_out\\tok\\t" >> ${og_id}.convergence_summary.tsv
        else
            echo -e "${og_id}\\tCSUBST\\t\\tfailed\\tsee_csubst.log" >> ${og_id}.convergence_summary.tsv
        fi
    else
        echo -e "${og_id}\\tCSUBST\\t\\tskipped\\tcsubst not installed" >> ${og_id}.convergence_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcoc: \$(command -v pcoc_detect >/dev/null && echo available || echo missing)
        csubst: \$(command -v csubst >/dev/null && echo available || echo missing)
    END_VERSIONS
    """
}
