process IQTREE_ASR {
    tag "${og_id}"
    label 'process_medium'
    publishDir "${params.outdir}/10_graft_phenotype/asr/${og_id}", mode: 'copy'

    input:
    tuple val(og_id), path(alignment)
    path shift_branches

    output:
    tuple val(og_id), path("${og_id}.asr_substitutions.tsv"), emit: substitutions
    path "${og_id}.treefile", emit: tree
    path "${og_id}.state", emit: state, optional: true
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    prep_og_alignment.py \\
        --alignment ${alignment} \\
        --out-fasta ${og_id}.species.faa \\
        --out-species-map ${og_id}.species_map.tsv

    NSEQ=\$(grep -c '^>' ${og_id}.species.faa || true)
    if [ "\$NSEQ" -lt 3 ]; then
        echo -e "Orthogroup\\tparent_node\\tchild_node\\tchild_tips\\tsite\\tparent_aa\\tchild_aa\\tparent_prob\\tchild_prob\\ton_phenotype_shift_branch\\tsubstitution" \\
            > ${og_id}.asr_substitutions.tsv
        echo "(${og_id}) skipped ASR: need >=3 species" >&2
        touch ${og_id}.treefile
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            iqtree: skipped
        END_VERSIONS
        exit 0
    fi

    IQTREE_BIN=\$(command -v iqtree2 || command -v iqtree)
    \$IQTREE_BIN \\
        -s ${og_id}.species.faa \\
        -m MFP \\
        -bb 1000 \\
        -asr \\
        -nt ${task.cpus} \\
        --prefix ${og_id} \\
        -redo

    STATE=\$(ls ${og_id}.state 2>/dev/null || true)
    TREE=${og_id}.treefile
    if [ -n "\$STATE" ] && [ -s "\$TREE" ]; then
        parse_iqtree_asr.py \\
            --state "\$STATE" \\
            --tree "\$TREE" \\
            --orthogroup ${og_id} \\
            --shift-branches ${shift_branches} \\
            --out-tsv ${og_id}.asr_substitutions.tsv \\
            --min-prob ${params.graft_asr_min_prob}
    else
        echo -e "Orthogroup\\tparent_node\\tchild_node\\tchild_tips\\tsite\\tparent_aa\\tchild_aa\\tparent_prob\\tchild_prob\\ton_phenotype_shift_branch\\tsubstitution" \\
            > ${og_id}.asr_substitutions.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(\$IQTREE_BIN --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
