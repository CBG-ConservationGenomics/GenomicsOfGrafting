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
    ch_codon_alns      // tuple(og_id, codon_aln) — optional empty

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
        // Queue channels can be read once. Split the trees when both consumers run.
        if (params.run_convergence || params.run_hyphy) {
            ch_tree_pairs = IQTREE_ASR.out.tree.map { t -> tuple(t.simpleName, t) }
        }
    } else {
        ch_asr = Channel.fromPath("${projectDir}/assets/empty.tsv").collect()
        ch_tree_pairs = Channel.empty()
    }

    if (params.run_convergence && params.run_hyphy && params.run_asr) {
        ch_tree_pairs.multiMap { og, tree ->
            convergence: tuple(og, tree)
            hyphy: tuple(og, tree)
        }.set { ch_trees }
        ch_tree_conv = ch_trees.convergence
        ch_tree_hyphy = ch_trees.hyphy
    } else if (params.run_convergence) {
        ch_tree_conv = ch_tree_pairs
        ch_tree_hyphy = Channel.empty()
    } else if (params.run_hyphy) {
        ch_tree_conv = Channel.empty()
        ch_tree_hyphy = ch_tree_pairs
    }

    if (params.run_convergence) {
        ch_conv_in = ch_og_alns.join(ch_tree_conv, by: 0, remainder: true)
            .map { og, aln, tree ->
                tree ? tuple(og, aln, tree) : null
            }
            .filter { it != null }

        RUN_CONVERGENCE(
            ch_conv_in,
            DERIVE_GRAFT_TRAITS.out.foreground.first()
        )
        ch_conv = RUN_CONVERGENCE.out.summary.map { _og, f -> f }.collect()
    } else {
        ch_conv = Channel.fromPath("${projectDir}/assets/empty.tsv").collect()
    }

    if (params.run_hyphy) {
        if (!params.run_asr) {
            error "--run_hyphy true requires --run_asr true so each orthogroup has a species tree"
        }
        HYPHY_SELECTION(
            ch_codon_alns.join(ch_tree_hyphy),
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
