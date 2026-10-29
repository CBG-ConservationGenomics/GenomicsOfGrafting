process HYPHY_SELECTION {
    tag "${og_id}"
    label 'process_high'
    publishDir "${params.outdir}/10_graft_phenotype/selection/${og_id}", mode: 'copy'

    input:
    tuple val(og_id), path(codon_alignment), path(tree)
    path foreground_tips

    output:
    tuple val(og_id), path("${og_id}.hyphy_summary.tsv"), emit: summary
    path "*.json", emit: json, optional: true
    path "${og_id}.fg.nwk", emit: tree, optional: true
    path "*.log", emit: logs, optional: true
    path "versions.yml", emit: versions

    when:
    params.run_hyphy

    script:
    // Shell locals are escaped (\$method). Nextflow would otherwise treat
    // them as its own variables and abort the task before HyPhy starts.
    """
    set -euo pipefail

    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' Orthogroup method pvalue key_stat status > ${og_id}.hyphy_summary.tsv

    prepare_hyphy_inputs.py prepare \\
        --codon ${codon_alignment} \\
        --tree ${tree} \\
        --foreground ${foreground_tips} \\
        --out-fasta ${og_id}.hyphy.fna \\
        --out-tree ${og_id}.fg.nwk \\
        --out-status ${og_id}.prep.tsv

    STATUS=\$(cut -f1 ${og_id}.prep.tsv)
    NOTE=\$(cut -f2 ${og_id}.prep.tsv)
    N_FG=\$(cut -f3 ${og_id}.prep.tsv)
    N_SEQ=\$(cut -f4 ${og_id}.prep.tsv)

    write_row() {
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$1" "\$2" "\$3" "\$4" "\$5" >> ${og_id}.hyphy_summary.tsv
    }

    if [ "\$STATUS" != "ok" ]; then
        write_row "${og_id}" HYPHY "" "\$NOTE" skipped
    elif ! command -v hyphy >/dev/null 2>&1; then
        write_row "${og_id}" HYPHY "" "hyphy not on PATH" skipped
    else
        BRANCH_ARGS=""
        if [ "\$N_FG" -gt 0 ]; then
            BRANCH_ARGS="--branches Foreground"
        fi

        # set -e is off for the method calls. One HyPhy failure is a table row,
        # not a failed Nextflow task. CPU is a HyPhy global and must precede
        # the method name. Stdin is closed so a missing option cannot prompt.
        set +e
        run_one() {
            local label="\$1"
            local method="\$2"
            local outjson="\$3"
            local extra="\$4"
            hyphy CPU=${task.cpus} "\$method" \\
                --alignment ${og_id}.hyphy.fna \\
                --tree ${og_id}.fg.nwk \\
                \$extra \\
                --output "\$outjson" \\
                < /dev/null > "\$method.log" 2>&1
            local rc=\$?
            if [ "\$rc" -eq 0 ] && [ -s "\$outjson" ]; then
                local pv
                pv=\$(prepare_hyphy_inputs.py pvalue "\$outjson" || true)
                write_row "${og_id}" "\$label" "\$pv" "" ok
            else
                local note
                note=\$(tail -n 1 "\$method.log" 2>/dev/null | sed 's/[[:space:]]\\+/ /g' | cut -c1-180 || true)
                write_row "${og_id}" "\$label" "" "\$note" failed
            fi
        }
        run_one BUSTED busted ${og_id}.BUSTED.json "\$BRANCH_ARGS"
        run_one aBSREL absrel ${og_id}.aBSREL.json "\$BRANCH_ARGS"
        run_one MEME meme ${og_id}.MEME.json "\$BRANCH_ARGS"
        if [ "\$N_FG" -gt 0 ] && [ "\$N_FG" -lt "\$N_SEQ" ]; then
            run_one RELAX relax ${og_id}.RELAX.json "--test Foreground"
        else
            write_row "${og_id}" RELAX "" "need both foreground and background branches" skipped
        fi
        set -e
    fi

    set +e
    ver=\$(hyphy --version </dev/null 2>/dev/null | head -n 1)
    set -e
    if [ -z "\$ver" ]; then
        ver="NA"
    fi
    {
        printf '"%s":\\n' "${task.process}"
        printf '    hyphy: %s\\n' "\$ver"
    } > versions.yml
    """
}
