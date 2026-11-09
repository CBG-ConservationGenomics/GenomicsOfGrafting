process RUN_CONVERGENCE {
    tag "${og_id}"
    label 'process_medium'
    publishDir "${params.outdir}/10_graft_phenotype/convergence/${og_id}", mode: 'copy'

    input:
    tuple val(og_id), path(alignment), path(tree)
    path foreground_tips
    path codon_dir

    output:
    tuple val(og_id), path("${og_id}.convergence_summary.tsv"), emit: summary
    path "${og_id}.pcoc_scenario.txt", optional: true
    path "pcoc_out", optional: true
    path "csubst_search", optional: true
    path "csubst_out", optional: true
    path "*.log", emit: logs, optional: true
    path "versions.yml", emit: versions

    when:
    params.run_convergence

    script:
    def pcoc_image = params.pcoc_image
    def csubst_image = params.csubst_image
    """
    set -euo pipefail

    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' Orthogroup method sites_or_stat status notes > ${og_id}.convergence_summary.tsv

    write_row() {
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$1" "\$2" "\$3" "\$4" "\$5" >> ${og_id}.convergence_summary.tsv
    }
    log_note() {
        tail -n 1 "\$1" 2>/dev/null | tr '\\t' ' ' | cut -c1-160 || true
    }

    prep_og_alignment.py \\
        --alignment ${alignment} \\
        --out-fasta ${og_id}.species.faa \\
        --out-species-map ${og_id}.species_map.tsv

    CODON_FAA=""
    if [ -d "${codon_dir}" ] && [ -s "${codon_dir}/${og_id}.codon.fna" ]; then
        prep_og_alignment.py \\
            --alignment "${codon_dir}/${og_id}.codon.fna" \\
            --out-fasta ${og_id}.species.codon.fna \\
            --out-species-map ${og_id}.codon_map.tsv
        CODON_FAA="${og_id}.species.codon.fna"
    fi

    # A host binary wins. Otherwise pull the image with Docker, or with
    # Singularity/Apptainer on clusters that do not expose the Docker socket.
    TOOL_ENGINE=""
    run_image() {
        local image="\$1"
        shift
        if command -v docker >/dev/null 2>&1; then
            local docker_ok=0
            if command -v timeout >/dev/null 2>&1; then
                timeout 20 docker info >/dev/null 2>&1 && docker_ok=1
            elif docker info >/dev/null 2>&1; then
                docker_ok=1
            fi
            if [ "\$docker_ok" -eq 1 ]; then
                TOOL_ENGINE="docker \$image"
                docker image inspect "\$image" >/dev/null 2>&1 || docker pull "\$image"
                docker run --rm \\
                    -e LOCAL_USER_ID="\$(id -u)" \\
                    -v "\$PWD:\$PWD" \\
                    -e CWD="\$PWD" \\
                    -w "\$PWD" \\
                    "\$image" "\$@"
                return
            fi
        fi
        local sbin=""
        if command -v singularity >/dev/null 2>&1; then
            sbin=singularity
        elif command -v apptainer >/dev/null 2>&1; then
            sbin=apptainer
        else
            echo "no docker daemon and no singularity/apptainer" >&2
            return 127
        fi
        TOOL_ENGINE="\$sbin \$image"
        "\$sbin" exec --bind "\$PWD" --pwd "\$PWD" "docker://\$image" "\$@"
    }

    if ! prep_convergence_inputs.py prepare \\
        --tree ${tree} \\
        --foreground ${foreground_tips} \\
        --protein-fasta ${og_id}.species.faa \\
        --out-tree ${og_id}.rooted.nwk \\
        --out-foreground ${og_id}.foreground.txt \\
        --out-status ${og_id}.prep.tsv \\
        \${CODON_FAA:+--codon-fasta "\$CODON_FAA"} \\
        --out-codon-tree ${og_id}.codon.nwk \\
        --out-codon-foreground ${og_id}.codon.foreground.txt \\
        > ${og_id}.prep.log 2>&1; then
        NOTE=\$(log_note ${og_id}.prep.log)
        write_row "${og_id}" PCOC "" "skipped" "\$NOTE"
        write_row "${og_id}" CSUBST "" "skipped" "\$NOTE"
    else
        STATUS=\$(awk -F '\\t' '\$1=="status"{print \$2}' ${og_id}.prep.tsv)
        NOTE=\$(awk -F '\\t' '\$1=="note"{print \$2}' ${og_id}.prep.tsv)
        N_FG=\$(awk -F '\\t' '\$1=="n_foreground"{print \$2}' ${og_id}.prep.tsv)
        N_CODON=\$(awk -F '\\t' '\$1=="n_codon"{print \$2}' ${og_id}.prep.tsv)
        N_CODON_FG=\$(awk -F '\\t' '\$1=="n_codon_foreground"{print \$2}' ${og_id}.prep.tsv)

        if [ "\$STATUS" != "ok" ]; then
            write_row "${og_id}" PCOC "" "skipped" "\$NOTE"
            write_row "${og_id}" CSUBST "" "skipped" "\$NOTE"
        else
            PCOC_ENGINE=""
            if command -v pcoc_det.py >/dev/null 2>&1 && command -v pcoc_num_tree.py >/dev/null 2>&1; then
                PCOC_ENGINE="host"
            fi
            if [ "\$N_FG" -lt 2 ]; then
                write_row "${og_id}" PCOC "" "skipped" "need at least two foreground species on this tree"
            elif [ "\$PCOC_ENGINE" = "host" ]; then
                if pcoc_num_tree.py -t ${og_id}.rooted.nwk -o ${og_id}.pcoc_num.pdf -n -u ${og_id}.pcoc_numbered.nwk > ${og_id}.pcoc.log 2>&1 \\
                    && prep_convergence_inputs.py scenario \\
                        --numbered ${og_id}.pcoc_numbered.nwk \\
                        --foreground ${og_id}.foreground.txt \\
                        --out-scenario ${og_id}.pcoc_scenario.txt >> ${og_id}.pcoc.log 2>&1 \\
                    && [ "\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt || true)" -ge 2 ] \\
                    && pcoc_det.py -t ${og_id}.rooted.nwk -aa ${og_id}.species.faa -o pcoc_out -m "\$(tr -d '[:space:]' < ${og_id}.pcoc_scenario.txt)" -f 0.8 >> ${og_id}.pcoc.log 2>&1; then
                    write_row "${og_id}" PCOC "n_events=\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt)" "ok" "host"
                elif [ -s ${og_id}.pcoc_scenario.txt ] && [ "\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt)" -lt 2 ]; then
                    write_row "${og_id}" PCOC "" "skipped" "foreground species did not map to two PCOC transitions"
                else
                    write_row "${og_id}" PCOC "" "failed" "host: \$(log_note ${og_id}.pcoc.log)"
                fi
            elif run_image "${pcoc_image}" pcoc_num_tree.py -t ${og_id}.rooted.nwk -o ${og_id}.pcoc_num.pdf -n -u ${og_id}.pcoc_numbered.nwk > ${og_id}.pcoc.log 2>&1 \\
                && prep_convergence_inputs.py scenario \\
                    --numbered ${og_id}.pcoc_numbered.nwk \\
                    --foreground ${og_id}.foreground.txt \\
                    --out-scenario ${og_id}.pcoc_scenario.txt >> ${og_id}.pcoc.log 2>&1 \\
                && [ "\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt || true)" -ge 2 ] \\
                && run_image "${pcoc_image}" pcoc_det.py -t ${og_id}.rooted.nwk -aa ${og_id}.species.faa -o pcoc_out -m "\$(tr -d '[:space:]' < ${og_id}.pcoc_scenario.txt)" -f 0.8 >> ${og_id}.pcoc.log 2>&1; then
                write_row "${og_id}" PCOC "n_events=\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt)" "ok" "\$TOOL_ENGINE"
            elif [ -s ${og_id}.pcoc_scenario.txt ] && [ "\$(awk -F/ 'NF{print NF; exit}' ${og_id}.pcoc_scenario.txt)" -lt 2 ]; then
                write_row "${og_id}" PCOC "" "skipped" "foreground species did not map to two PCOC transitions"
            else
                write_row "${og_id}" PCOC "" "failed" "\${TOOL_ENGINE:-image}: \$(log_note ${og_id}.pcoc.log)"
            fi

            CSUBST_ENGINE=""
            if command -v csubst >/dev/null 2>&1; then
                CSUBST_ENGINE="host"
            fi
            if [ "\$N_CODON" -lt 3 ]; then
                write_row "${og_id}" CSUBST "" "skipped" "no codon alignment for this orthogroup"
            elif [ "\$N_CODON_FG" -lt 1 ]; then
                write_row "${og_id}" CSUBST "" "skipped" "no foreground species in the codon alignment"
            elif [ "\$CSUBST_ENGINE" = "host" ]; then
                if csubst search \\
                    --alignment_file ${og_id}.species.codon.fna \\
                    --rooted_tree_file ${og_id}.codon.nwk \\
                    --foreground ${og_id}.codon.foreground.txt \\
                    --outdir csubst_out \\
                    > ${og_id}.csubst.log 2>&1; then
                    write_row "${og_id}" CSUBST "see_csubst_out" "ok" "host"
                else
                    write_row "${og_id}" CSUBST "" "failed" "host: \$(log_note ${og_id}.csubst.log)"
                fi
            elif run_image "${csubst_image}" csubst search \\
                --alignment_file ${og_id}.species.codon.fna \\
                --rooted_tree_file ${og_id}.codon.nwk \\
                --foreground ${og_id}.codon.foreground.txt \\
                --outdir csubst_out \\
                > ${og_id}.csubst.log 2>&1; then
                write_row "${og_id}" CSUBST "see_csubst_out" "ok" "\$TOOL_ENGINE"
            else
                write_row "${og_id}" CSUBST "" "failed" "\${TOOL_ENGINE:-image}: \$(log_note ${og_id}.csubst.log)"
            fi
        fi
    fi

    {
        printf '"%s":\\n' "${task.process}"
        printf '    pcoc: %s\\n' "${pcoc_image}"
        printf '    csubst: %s\\n' "${csubst_image}"
    } > versions.yml
    """
}
