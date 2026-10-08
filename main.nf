#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * GraftingGeneFamilies
 * Download plant proteomes → QC/rename → BUSCO gate → OrthoFinder →
 * select GOI orthogroups → annotate selected members →
 * protein % identity → polished tables
 */

include { DOWNLOAD_PROTEOME           } from './modules/local/download_proteome'
include { QC_PROTEOME                 } from './modules/local/qc_proteome'
include { QC_CDS                      } from './modules/local/qc_cds'
include { BUILD_CODON_ALIGNMENTS      } from './modules/local/build_codon_alignments'
include { RUN_BUSCO                   } from './modules/local/run_busco'
include { FILTER_BUSCO                } from './modules/local/filter_busco'
include { COLLECT_BUSCO_SUMMARY       } from './modules/local/collect_busco_summary'
include { RUN_ORTHOFINDER             } from './modules/local/run_orthofinder'
include { SELECT_ORTHOGROUPS          } from './modules/local/select_orthogroups'
include { EXTRACT_SELECTED_PROTEINS   } from './modules/local/extract_selected_proteins'
include { ANNOTATE_EGGNOG             } from './modules/local/annotate_eggnog'
include { ANNOTATE_INTERPROSCAN       } from './modules/local/annotate_interproscan'
include { ANNOTATE_SIGNALP            } from './modules/local/annotate_signalp'
include { ANNOTATE_DEEPTMHMM          } from './modules/local/annotate_deeptmhmm'
include { CALC_PERCENT_IDENTITY       } from './modules/local/calc_percent_identity'
include { MAKE_REPORT                 } from './modules/local/make_report'
include { PREPARE_ALPHAFOLD_INPUTS    } from './modules/local/prepare_alphafold_inputs'
include { GRAFT_PHENOTYPE             } from './subworkflows/local/graft_phenotype'

def helpMessage() {
    log.info """
    GraftingGeneFamilies v${workflow.manifest.version}

    Usage:
      nextflow run main.nf -profile hpc,conda \\
        --genes_of_interest data/genes_of_interest.txt \\
        --outdir results

    Key parameters:
      --species_yaml          Species list (default: conf/species.yaml)
      --genes_of_interest     One gene ID per line (e.g. Solyc01g000570)
      --outdir                Output directory
      --busco_min_complete    Min Complete BUSCO % to keep a species (default: 90)
      --busco_lineage         Default lineage (default: embryophyta_odb10)
      --busco_download_path   Shared BUSCO download/cache directory on the server
      --reference_species_id  %ID reference: auto (matched GOI) or Sly/Ath/Nbe/Cpe
      --run_eggnog / --run_interproscan / --run_signalp / --run_deeptmhmm
      --eggnog_data_dir       Path to eggNOG DB on the server
      --interproscan_dir      Path to InterProScan installation
    """.stripIndent()
}

def loadSpecies(String yamlPath) {
    def text = new File(yamlPath).text
    def species = []
    def current = null
    text.split('\n').each { line ->
        def t = line.replaceAll(/#.*/, '').stripIndent().trim()
        if (t.startsWith('- id:')) {
            if (current) species << current
            current = [:]
            current.id = t.replaceFirst(/^- id:\s*/, '').replaceAll(/^["']|["']$/, '')
        } else if (current != null && t.contains(':')) {
            def parts = t.split(':', 2)
            def key = parts[0].trim()
            def val = parts[1].trim().replaceAll(/^["']|["']$/, '')
            if (val == 'null' || val == '~' || val == '') {
                current[key] = null
            } else {
                current[key] = val
            }
        }
    }
    if (current) species << current
    return species
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    def species_list = loadSpecies(params.species_yaml)
    if (!species_list) {
        error "No species parsed from ${params.species_yaml}"
    }

    Channel
        .fromList(species_list)
        .map { s ->
            [
                id            : s.id,
                genus         : s.genus,
                species       : s.species,
                source        : s.source,
                accession     : s.accession,
                fasta_url     : s.fasta_url,
                cds_url       : s.cds_url,
                gff_url       : s.gff_url,
                fasta_path    : s.fasta_path,
                cds_path      : s.cds_path,
                gff_path      : s.gff_path,
                busco_lineage : s.busco_lineage ?: params.busco_lineage,
                notes         : s.notes
            ]
        }
        .set { ch_species }

    // Download proteins + CDS → QC/rename
    DOWNLOAD_PROTEOME(ch_species)

    DOWNLOAD_PROTEOME.out.fasta
        .join(DOWNLOAD_PROTEOME.out.gff, by: 0)
        .set { ch_raw }

    QC_PROTEOME(ch_raw)

    // Align CDS to retained proteins (same IDs); trim ORFs when transcripts include UTRs
    DOWNLOAD_PROTEOME.out.cds
        .join(QC_PROTEOME.out.fasta, by: 0)
        .join(QC_PROTEOME.out.id_map, by: 0)
        .set { ch_cds_qc_in }

    QC_CDS(ch_cds_qc_in)

    // -------------------------------------------------------------------------
    // BUSCO completeness gate — failing species are dropped before OrthoFinder
    // -------------------------------------------------------------------------
    if (params.run_busco) {
        RUN_BUSCO(QC_PROTEOME.out.fasta)

        QC_PROTEOME.out.fasta
            .join(QC_PROTEOME.out.id_map, by: 0)
            .join(RUN_BUSCO.out.metrics, by: 0)
            .join(RUN_BUSCO.out.status, by: 0)
            .map { meta, fasta, id_map, metrics, status ->
                [meta, fasta, id_map, metrics, status]
            }
            .set { ch_busco_gate_in }

        FILTER_BUSCO(ch_busco_gate_in)

        COLLECT_BUSCO_SUMMARY(
            RUN_BUSCO.out.metrics.map { _meta, m -> m }.collect()
        )

        FILTER_BUSCO.out.pass_fasta
            .toList()
            .map { rows ->
                def n = rows.size()
                def min_n = params.busco_min_species as int
                if (n < min_n) {
                    error "BUSCO gate: only ${n} species passed (>= ${params.busco_min_complete}% Complete). Need >= ${min_n}. See ${params.outdir}/02b_busco/"
                }
                log.info "BUSCO gate: ${n} species passed (>= ${params.busco_min_complete}% Complete)"
                rows
            }
            .flatMap { it }
            .set { ch_pass_fasta }

        FILTER_BUSCO.out.pass_id_map
            .map { _meta, m -> m }
            .collect()
            .set { ch_id_maps }

        ch_pass_fasta
            .map { _meta, fa -> fa }
            .collect()
            .set { ch_proteomes }
    } else {
        log.warn "BUSCO gate disabled (--run_busco false); all QC'd proteomes proceed"
        QC_PROTEOME.out.fasta
            .map { _meta, fa -> fa }
            .collect()
            .set { ch_proteomes }

        QC_PROTEOME.out.id_map
            .map { _meta, m -> m }
            .collect()
            .set { ch_id_maps }
    }

    // OrthoFinder (MSA mode so alignments are available for %ID)
    RUN_ORTHOFINDER(ch_proteomes)

    // Select orthogroups containing genes of interest
    ch_goi = Channel.fromPath(params.genes_of_interest, checkIfExists: true)

    SELECT_ORTHOGROUPS(
        RUN_ORTHOFINDER.out.orthogroups,
        ch_goi,
        ch_id_maps
    )

    // Protein-MSA-guided codon alignments (for HyPhy / selection)
    if (params.build_codon_alignments) {
        BUILD_CODON_ALIGNMENTS(
            SELECT_ORTHOGROUPS.out.list,
            RUN_ORTHOFINDER.out.alignments,
            QC_CDS.out.cds.map { _meta, cds -> cds }.collect()
        )
        ch_codon_dir = BUILD_CODON_ALIGNMENTS.out.dir
    } else {
        ch_codon_dir = Channel.empty()
    }

    // Extract only selected OG proteins for annotation
    EXTRACT_SELECTED_PROTEINS(
        SELECT_ORTHOGROUPS.out.tsv,
        ch_proteomes
    )

    ch_ann_fasta = EXTRACT_SELECTED_PROTEINS.out.fasta

    if (params.run_eggnog) {
        ANNOTATE_EGGNOG(ch_ann_fasta)
        ch_eggnog = ANNOTATE_EGGNOG.out.annotations
    } else {
        ch_eggnog = Channel.fromPath("${projectDir}/assets/empty.tsv")
    }

    if (params.run_interproscan) {
        ANNOTATE_INTERPROSCAN(ch_ann_fasta)
        ch_ips = ANNOTATE_INTERPROSCAN.out.tsv
    } else {
        ch_ips = Channel.fromPath("${projectDir}/assets/empty.tsv")
    }

    if (params.run_signalp) {
        ANNOTATE_SIGNALP(ch_ann_fasta)
        ch_signalp = ANNOTATE_SIGNALP.out.tsv
    } else {
        ch_signalp = Channel.fromPath("${projectDir}/assets/empty.tsv")
    }

    if (params.run_deeptmhmm) {
        ANNOTATE_DEEPTMHMM(ch_ann_fasta)
        ch_tm = ANNOTATE_DEEPTMHMM.out.tsv
    } else {
        ch_tm = Channel.fromPath("${projectDir}/assets/empty.tsv")
    }

    CALC_PERCENT_IDENTITY(
        SELECT_ORTHOGROUPS.out.list,
        SELECT_ORTHOGROUPS.out.tsv,
        RUN_ORTHOFINDER.out.alignments
    )

    MAKE_REPORT(
        SELECT_ORTHOGROUPS.out.tsv,
        CALC_PERCENT_IDENTITY.out.summary,
        CALC_PERCENT_IDENTITY.out.pairwise,
        ch_eggnog,
        ch_ips,
        ch_signalp,
        ch_tm
    )

    if (params.prepare_alphafold) {
        PREPARE_ALPHAFOLD_INPUTS(
            SELECT_ORTHOGROUPS.out.tsv,
            CALC_PERCENT_IDENTITY.out.summary,
            CALC_PERCENT_IDENTITY.out.pairwise,
            EXTRACT_SELECTED_PROTEINS.out.fasta
        )
    }

    // Graft outcome matrix → traits, ASR, divergence tests, optional HyPhy/convergence
    if (params.run_graft_phenotype) {
        ch_og_alns = SELECT_ORTHOGROUPS.out.list
            .splitText()
            .map { it.trim() }
            .filter { it }
            .combine(RUN_ORTHOFINDER.out.alignments)
            .map { og, alndir ->
                def cand = ['.fa', '.faa', '.fasta', '.aln']
                    .collect { ext -> file("${alndir}/${og}${ext}") }
                    .find { it.exists() }
                cand ? tuple(og, cand) : null
            }
            .filter { it != null }

        ch_codon = Channel.empty()
        if (params.run_hyphy) {
            if (params.build_codon_alignments) {
                ch_codon = BUILD_CODON_ALIGNMENTS.out.dir
                    .flatMap { dir ->
                        def list = file("${dir}/*.codon.fna")
                        list instanceof List ? list : [list]
                    }
                    .filter { it && it.exists() }
                    .map { f ->
                        def og = f.name.replaceAll(/\.codon\.fna$/, '')
                        tuple(og, f, file("${params.outdir}/10_graft_phenotype/asr/${og}/${og}.treefile"))
                    }
            } else if (params.codon_alignments_dir) {
                ch_codon = Channel
                    .fromPath("${params.codon_alignments_dir}/*.{fa,fasta,fna}")
                    .map { f -> tuple(f.simpleName.replaceAll(/\.codon$/, ''), f, file('NO_TREE')) }
            } else {
                error "--run_hyphy true requires --build_codon_alignments true or --codon_alignments_dir"
            }
        }

        GRAFT_PHENOTYPE(
            Channel.fromPath(params.graft_matrix, checkIfExists: true),
            RUN_ORTHOFINDER.out.species_tree,
            SELECT_ORTHOGROUPS.out.list,
            RUN_ORTHOFINDER.out.alignments,
            ch_og_alns,
            ch_codon
        )
    }
}

workflow.onComplete {
    log.info """
    Pipeline complete.
      Success : ${workflow.success}
      Outdir  : ${params.outdir}
      BUSCO   : ${params.outdir}/02b_busco/
      Tables  : ${params.outdir}/07_tables/
    """.stripIndent()
}
