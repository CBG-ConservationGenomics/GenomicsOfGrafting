process DOWNLOAD_PROTEOME {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/01_raw_proteomes", mode: 'copy'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.proteins.faa"), emit: fasta
    tuple val(meta), path("${meta.id}.cds.fna"), emit: cds
    tuple val(meta), path("${meta.id}.gff3"), emit: gff
    path "versions.yml", emit: versions

    script:
    def accession  = meta.accession ?: ''
    def fasta_url  = meta.fasta_url ?: ''
    def cds_url    = meta.cds_url ?: ''
    def gff_url    = meta.gff_url ?: ''
    def fasta_path = meta.fasta_path ?: ''
    def cds_path   = meta.cds_path ?: ''
    def gff_path   = meta.gff_path ?: ''
    """
    set -euo pipefail

    stage_file() {
        local src="\$1"
        local dest="\$2"
        if [ ! -s "\$src" ]; then
            echo "ERROR: missing or empty source: \$src" >&2
            exit 1
        fi
        local magic
        magic=\$(od -An -tx1 -N2 "\$src" | tr -d ' \\n')
        if [ "\$magic" = "1f8b" ]; then
            gzip -dc "\$src" > "\$dest"
        elif [ "\$magic" = "425a" ]; then
            bzip2 -dc "\$src" > "\$dest"
        else
            cp "\$src" "\$dest"
        fi
    }

    validate_fasta() {
        local f="\$1"
        local label="\$2"
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$label is empty" >&2
            exit 1
        fi
        local magic
        magic=\$(od -An -tx1 -N2 "\$f" | tr -d ' \\n')
        if [ "\$magic" = "1f8b" ]; then
            echo "ERROR: \$label is still gzip-compressed" >&2
            exit 1
        fi
        if head -n1 "\$f" | grep -qiE '<(!DOCTYPE|html)|Not Found'; then
            echo "ERROR: \$label looks like HTML" >&2
            head -n5 "\$f" >&2
            exit 1
        fi
        if ! head -n1 "\$f" | grep -q '^>'; then
            echo "ERROR: \$label does not start with '>'" >&2
            head -n5 "\$f" >&2
            exit 1
        fi
    }

    if [ "${meta.source}" = "ncbi" ]; then
        datasets download genome accession ${accession} \\
            --include protein,gff3,cds \\
            --filename ${meta.id}_ncbi.zip
        unzip -q ${meta.id}_ncbi.zip -d ${meta.id}_ncbi
        fa=\$(find ${meta.id}_ncbi -name 'protein.faa' | head -n1)
        cds=\$(find ${meta.id}_ncbi \\( -name 'cds_from_genomic.fna' -o -name '*cds*.fna' -o -name '*cds*.fa' \\) | head -n1)
        gff=\$(find ${meta.id}_ncbi \\( -name 'genomic.gff' -o -name '*.gff3' -o -name '*.gff' \\) | head -n1)
        cp "\$fa" ${meta.id}.proteins.faa
        if [ -n "\$cds" ]; then
            cp "\$cds" ${meta.id}.cds.fna
        else
            echo "WARNING: no CDS in NCBI package for ${meta.id}; creating empty placeholder" >&2
            touch ${meta.id}.cds.fna
        fi
        if [ -n "\$gff" ]; then
            cp "\$gff" ${meta.id}.gff3
        else
            touch ${meta.id}.gff3
        fi
    elif [ "${meta.source}" = "url" ]; then
        curl -L --fail --retry 3 -o ${meta.id}.proteins.download "${fasta_url}"
        stage_file ${meta.id}.proteins.download ${meta.id}.proteins.faa
        if [ -n "${cds_url}" ] && [ "${cds_url}" != "null" ]; then
            curl -L --fail --retry 3 -o ${meta.id}.cds.download "${cds_url}"
            stage_file ${meta.id}.cds.download ${meta.id}.cds.fna
        else
            echo "WARNING: no cds_url for ${meta.id}" >&2
            touch ${meta.id}.cds.fna
        fi
        if [ -n "${gff_url}" ] && [ "${gff_url}" != "null" ]; then
            curl -L --fail --retry 3 -o ${meta.id}.gff.download "${gff_url}" \\
                && stage_file ${meta.id}.gff.download ${meta.id}.gff3 \\
                || touch ${meta.id}.gff3
        else
            touch ${meta.id}.gff3
        fi
    elif [ "${meta.source}" = "local" ]; then
        stage_file "${fasta_path}" ${meta.id}.proteins.faa
        if [ -n "${cds_path}" ] && [ -f "${cds_path}" ]; then
            stage_file "${cds_path}" ${meta.id}.cds.fna
        else
            touch ${meta.id}.cds.fna
        fi
        if [ -n "${gff_path}" ] && [ -f "${gff_path}" ]; then
            stage_file "${gff_path}" ${meta.id}.gff3
        else
            touch ${meta.id}.gff3
        fi
    else
        echo "Unknown source: ${meta.source}" >&2
        exit 1
    fi

    validate_fasta ${meta.id}.proteins.faa proteins
    # CDS may be empty for some species; only validate if non-empty
    if [ -s ${meta.id}.cds.fna ]; then
        validate_fasta ${meta.id}.cds.fna cds
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datasets: \$(datasets version 2>/dev/null | head -n1 || echo 'NA')
        curl: \$(curl --version | head -n1)
    END_VERSIONS
    """
}
