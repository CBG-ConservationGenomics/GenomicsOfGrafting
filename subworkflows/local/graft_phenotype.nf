include { DERIVE_GRAFT_TRAITS         } from '../../modules/local/derive_graft_traits'
include { GRAFT_PAIR_DIVERGENCE       } from '../../modules/local/graft_pair_divergence'
include { TEST_DIVERGENCE_PHENOTYPE   } from '../../modules/local/test_divergence_phenotype'
include { IQTREE_ASR                  } from '../../modules/local/iqtree_asr'
include { HYPHY_SELECTION             } from '../../modules/local/hyphy_selection'
include { RUN_CONVERGENCE             } from '../../modules/local/run_convergence'
include { RANK_GRAFT_CANDIDATES       } from '../../modules/local/rank_graft_candidates'

workflow GRAFT_PHENOTYPE {
    take:
    ch_matrix          // path
    ch_species_tree    // path
    ch_og_list         // path
    ch_alignments_dir  // path
    ch_og_alns         // tuple(og_id, alignment) — optional empty
    ch_codon_alns      // tuple(og_id, codon_aln, tree) — optional empty

    main:
    DERIVE_GRAFT_TRAITS(ch_matrix, ch_species_tree)

    GRAFT_PAIR_DIVERGENCE(
        DERIVE_GRAFT_TRAITS.out.pairs,
        ch_og_list,
        ch_alignments_dir
    )

    TEST_DIVERGENCE_PHENOTYPE(
        GRAFT_PAIR_DIVERGENCE.out.divergence,
        DERIVE_GRAFT_TRAITS.out.sp_distance
    )

    if (params.run_asr) {
        IQTREE_ASR(
            ch_og_alns,
            DERIVE_GRAFT_TRAITS.out.shift_branches.first()
        )
        ch_asr = IQTREE_ASR.out.substitutions.map { _og, f -> f }.collect()
        ch_trees = IQTREE_ASR.out.tree
    } else {
        ch_asr = Channel.fromPath("${projectDir}/assets/empty.tsv").collect()
        ch_trees = Channel.empty()
    }

    if (params.run_convergence) {
        ch_conv_in = ch_og_alns.join(
            ch_trees.map { t ->
                def og = t.simpleName
                tuple(og, t)
            },
            by: 0,
            remainder: true
        ).map { og, aln, tree ->
            tree ? tuple(og, aln, tree) : null
        }.filter { it != null }

        RUN_CONVERGENCE(
            ch_conv_in,
            DERIVE_GRAFT_TRAITS.out.foreground.first()
        )
        ch_conv = RUN_CONVERGENCE.out.summary.map { _og, f -> f }.collect()
    } else {
        ch_conv = Channel.fromPath("${projectDir}/assets/empty.tsv").collect()
    }

    if (params.run_hyphy) {
        HYPHY_SELECTION(
            ch_codon_alns,
            DERIVE_GRAFT_TRAITS.out.foreground.first()
        )
        ch_sel = HYPHY_SELECTION.out.summary.map { _og, f -> f }.collect()
    } else {
        ch_sel = Channel.fromPath("${projectDir}/assets/empty.tsv").collect()
    }

    RANK_GRAFT_CANDIDATES(
        ch_asr,
        TEST_DIVERGENCE_PHENOTYPE.out.tests,
        ch_sel,
        ch_conv
    )

    emit:
    pairs     = DERIVE_GRAFT_TRAITS.out.pairs
    lineage   = DERIVE_GRAFT_TRAITS.out.lineage
    ranking   = RANK_GRAFT_CANDIDATES.out.ranking
    div_tests = TEST_DIVERGENCE_PHENOTYPE.out.tests
}
