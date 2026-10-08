#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Standalone graft-phenotype workflow. See docs/graft_phenotype.md
 */

include { GRAFT_PHENOTYPE } from '../subworkflows/local/graft_phenotype'

params.graft_matrix           = null
params.species_tree           = null
params.og_list                = null
params.alignments_dir         = null
params.codon_alignments_dir   = null
params.outdir                 = 'results'

if (!params.graft_matrix) error "Set --graft_matrix"
if (!params.species_tree) error "Set --species_tree"
if (!params.og_list) error "Set --og_list"
if (!params.alignments_dir) error "Set --alignments_dir"

workflow {
    ch_matrix = Channel.fromPath(params.graft_matrix, checkIfExists: true)
    ch_stree  = Channel.fromPath(params.species_tree, checkIfExists: true)
    ch_oglist = Channel.fromPath(params.og_list, checkIfExists: true)
    ch_aln    = Channel.fromPath(params.alignments_dir, checkIfExists: true)

    ch_og_alns = Channel
        .fromPath(params.og_list, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { og ->
            def dir = file(params.alignments_dir)
            def cand = ['.fa', '.faa', '.fasta', '.aln']
                .collect { ext -> file("${dir}/${og}${ext}") }
                .find { it.exists() }
            cand ? tuple(og, cand) : null
        }
        .filter { it != null }

    if (params.run_hyphy) {
        if (!params.codon_alignments_dir) {
            error "--run_hyphy requires --codon_alignments_dir"
        }
        ch_codon = Channel
            .fromPath("${params.codon_alignments_dir}/*.{fa,fasta,fna}")
            .map { f ->
                def og = f.simpleName
                def tree = file("${params.outdir}/10_graft_phenotype/asr/${og}/${og}.treefile")
                tuple(og, f, tree)
            }
    } else {
        ch_codon = Channel.empty()
    }

    GRAFT_PHENOTYPE(
        ch_matrix,
        ch_stree,
        ch_oglist,
        ch_aln,
        ch_og_alns,
        ch_codon
    )
}

workflow.onComplete {
    log.info "Graft phenotype done → ${params.outdir}/10_graft_phenotype/"
}
