process QC_PROTEOME {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/02_clean_proteomes", mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(gff)

    output:
    tuple val(meta), path("${meta.id}.clean.faa"), emit: fasta
    tuple val(meta), path("${meta.id}.id_map.tsv"), emit: id_map
    path "versions.yml", emit: versions

    script:
    def isoform_flag = params.keep_longest_isoform ? '--keep-longest-isoform' : '--no-keep-longest-isoform'
    def hap_flag     = params.drop_alt_haplotigs ? '--drop-alt-haplotigs' : '--no-drop-alt-haplotigs'
    """
    set -euo pipefail

    GFF_ARG=""
    if [ -s "${gff}" ]; then
        GFF_ARG="--gff ${gff}"
    fi

    qc_and_rename_proteome.py \\
        --fasta ${fasta} \\
        --species-id ${meta.id} \\
        --genus ${meta.genus} \\
        --species ${meta.species} \\
        --out-fasta ${meta.id}.clean.faa \\
        --out-map ${meta.id}.id_map.tsv \\
        \$GFF_ARG \\
        ${isoform_flag} \\
        ${hap_flag} \\
        --haplotig-patterns '${params.haplotig_patterns}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
