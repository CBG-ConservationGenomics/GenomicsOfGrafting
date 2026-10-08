process QC_CDS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/02_clean_cds", mode: 'copy'

    input:
    tuple val(meta), path(cds), path(proteins), path(id_map)

    output:
    tuple val(meta), path("${meta.id}.clean.cds.fna"), emit: cds, optional: true
    path "${meta.id}.cds_qc_report.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    if [ ! -s "${cds}" ]; then
        echo -e "n_cds_input\\tn_proteins\\tn_cds_kept\\tn_exact_or_ok\\tn_trimmed_orf\\tn_orf_fail\\tn_no_protein_match\\tn_duplicate\\tn_proteins_missing_cds" \\
            > ${meta.id}.cds_qc_report.tsv
        echo -e "0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0" >> ${meta.id}.cds_qc_report.tsv
        echo "WARNING: empty CDS for ${meta.id}; skipping" >&2
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version | sed 's/Python //')
        END_VERSIONS
        exit 0
    fi

    qc_and_rename_cds.py \\
        --cds ${cds} \\
        --proteins ${proteins} \\
        --id-map ${id_map} \\
        --out-cds ${meta.id}.clean.cds.fna \\
        --out-report ${meta.id}.cds_qc_report.tsv \\
        --trim-to-protein

    if [ ! -s ${meta.id}.clean.cds.fna ]; then
        echo "WARNING: no CDS retained for ${meta.id}" >&2
        rm -f ${meta.id}.clean.cds.fna
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
