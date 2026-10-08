process FILTER_BUSCO {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/02b_busco/gate", mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(id_map), path(metrics), path(status)

    output:
    tuple val(meta), path("${meta.id}.clean.faa"), emit: pass_fasta, optional: true
    tuple val(meta), path("${meta.id}.id_map.tsv"), emit: pass_id_map, optional: true
    path "${meta.id}.busco_gate.log", emit: log
    path "${meta.id}.busco_metrics.tsv", emit: metrics

    script:
    """
    set -euo pipefail

    STATUS=\$(tr -d '[:space:]' < ${status})
    COMPLETE=\$(awk -F'\\t' 'NR==2 {print \$3}' ${metrics})
    THRESH=${params.busco_min_complete}

    cp ${metrics} ${meta.id}.busco_metrics.tsv

    {
      echo "species_id=${meta.id}"
      echo "busco_complete_pct=\${COMPLETE}"
      echo "threshold_pct=\${THRESH}"
      echo "lineage=${meta.busco_lineage ?: params.busco_lineage}"
      echo "status=\${STATUS}"
    } > ${meta.id}.busco_gate.log

    if [ "\$STATUS" = "pass" ]; then
        cp ${fasta} ${meta.id}.clean.faa
        cp ${id_map} ${meta.id}.id_map.tsv
        echo "GATE=PASS — ${meta.id} proceeds to OrthoFinder" >> ${meta.id}.busco_gate.log
    else
        echo "GATE=FAIL — ${meta.id} EXCLUDED from OrthoFinder and downstream analyses" >> ${meta.id}.busco_gate.log
        echo "Reason: complete BUSCOs \${COMPLETE}% < \${THRESH}% required" >> ${meta.id}.busco_gate.log
    fi
    """
}
