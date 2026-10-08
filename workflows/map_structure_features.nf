#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Map functional annotations onto AlphaFold structures + optional:
 *   - MSA transfer of curated active-site / pore residues
 *   - FreeSASA / Biopython surface for defense proteins
 *   - Multimer interface PAE / contacts
 *
 *   nextflow run workflows/map_structure_features.nf -profile hpc,conda \
 *     --proteins results/05_selected_orthogroups/selected_og_proteins.faa \
 *     --signalp results/04_annotation/signalp/signalp.tsv \
 *     --interproscan results/04_annotation/interproscan/interproscan.tsv \
 *     --deeptmhmm results/04_annotation/deeptmhmm/deeptmhmm.tsv \
 *     --eggnog results/04_annotation/eggnog/eggnog.emapper.annotations \
 *     --models_dir alphafold_results/models \
 *     --substitutions results/09_structure_comparison/substitutions.tsv \
 *     --alignments_dir results/03_orthofinder/OrthoFinder/.../MultipleSequenceAlignments \
 *     --selected_ogs results/05_selected_orthogroups/selected_orthogroups.tsv \
 *     --multimer_models_dir alphafold_multimer/models \
 *     --multimer_manifest results/08_alphafold_inputs/multimer/manifest.tsv \
 *     --outdir results
 */

include { TRANSFER_CURATED_SITES      } from '../modules/local/transfer_curated_sites'
include { EXTRACT_SEQUENCE_FEATURES   } from '../modules/local/extract_sequence_features'
include { MAP_FEATURES_TO_STRUCTURE   } from '../modules/local/map_features_to_structure'
include { COMPUTE_SURFACE_RSA         } from '../modules/local/compute_surface_rsa'
include { PARSE_MULTIMER_INTERFACES   } from '../modules/local/parse_multimer_interfaces'

params.proteins               = null
params.signalp                = "${projectDir}/assets/empty.tsv"
params.interproscan           = "${projectDir}/assets/empty.tsv"
params.deeptmhmm              = "${projectDir}/assets/empty.tsv"
params.eggnog                 = "${projectDir}/assets/empty.tsv"
params.rules                  = "${projectDir}/conf/functional_feature_rules.yaml"
params.models_dir             = null
params.substitutions          = "${projectDir}/assets/empty.tsv"
params.alignments_dir         = null
params.selected_ogs           = "${projectDir}/assets/empty.tsv"
params.multimer_models_dir    = null
params.multimer_manifest      = null
params.outdir                 = 'results'
params.af_min_plddt           = 50.0
params.rsa_surface_cutoff     = 0.25
params.interface_distance_cutoff = 5.0
params.interface_pae_cutoff   = 15.0
params.run_surface            = true
params.run_multimer_interfaces = true

if (!params.proteins) error "Set --proteins"
if (!params.models_dir) error "Set --models_dir (AlphaFold/ColabFold models)"

workflow {
    ch_empty = Channel.fromPath("${projectDir}/assets/empty.tsv")

    // 1) Transfer curated sites via OrthoFinder MSA (optional)
    if (params.alignments_dir) {
        TRANSFER_CURATED_SITES(
            Channel.fromPath(params.rules, checkIfExists: true),
            Channel.fromPath(params.alignments_dir, checkIfExists: true),
            Channel.fromPath(params.selected_ogs)
        )
        ch_xfer = TRANSFER_CURATED_SITES.out.transferred
    } else {
        ch_xfer = ch_empty
    }

    // 2) Sequence features (+ MSA-transferred curated sites when available)
    EXTRACT_SEQUENCE_FEATURES(
        Channel.fromPath(params.proteins, checkIfExists: true),
        Channel.fromPath(params.signalp),
        Channel.fromPath(params.interproscan),
        Channel.fromPath(params.deeptmhmm),
        Channel.fromPath(params.eggnog),
        Channel.fromPath(params.rules, checkIfExists: true),
        ch_xfer
    )

    MAP_FEATURES_TO_STRUCTURE(
        EXTRACT_SEQUENCE_FEATURES.out.features,
        Channel.fromPath(params.models_dir, checkIfExists: true),
        Channel.fromPath(params.substitutions)
    )

    // 3) Defense surface RSA
    if (params.run_surface) {
        COMPUTE_SURFACE_RSA(
            Channel.fromPath(params.models_dir, checkIfExists: true),
            EXTRACT_SEQUENCE_FEATURES.out.classifications,
            Channel.fromPath(params.substitutions)
        )
    }

    // 4) Multimer interface PAE
    if (params.run_multimer_interfaces && params.multimer_models_dir && params.multimer_manifest) {
        PARSE_MULTIMER_INTERFACES(
            Channel.fromPath(params.multimer_models_dir, checkIfExists: true),
            Channel.fromPath(params.multimer_manifest, checkIfExists: true)
        )
    }
}

workflow.onComplete {
    log.info "Feature mapping done → ${params.outdir}/09_structure_comparison/"
}
