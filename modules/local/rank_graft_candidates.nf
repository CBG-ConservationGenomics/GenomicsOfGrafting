process RANK_GRAFT_CANDIDATES {
    label 'process_low'
    publishDir "${params.outdir}/10_graft_phenotype", mode: 'copy'

    input:
    path asr_dir
    path divergence_tests
    path selection_dir
    path convergence_dir

    output:
    path "graft_candidate_ranking.tsv", emit: ranking
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    mkdir -p asr_in selection_in convergence_in
    # Flatten collected inputs (may be files or directories)
    cp -r ${asr_dir} asr_in/ 2>/dev/null || true
    find . -name '*.asr_substitutions.tsv' -exec cp {} asr_in/ \\; 2>/dev/null || true

    SEL_ARG=""
    CONV_ARG=""
    if [ -e "${selection_dir}" ]; then
        mkdir -p selection_in
        cp -r ${selection_dir} selection_in/ 2>/dev/null || true
        find . -name '*.hyphy_summary.tsv' -exec cp {} selection_in/ \\; 2>/dev/null || true
        SEL_ARG="--selection-dir selection_in"
    fi
    if [ -e "${convergence_dir}" ]; then
        mkdir -p convergence_in
        cp -r ${convergence_dir} convergence_in/ 2>/dev/null || true
        find . -name '*.convergence_summary.tsv' -exec cp {} convergence_in/ \\; 2>/dev/null || true
        CONV_ARG="--convergence-dir convergence_in"
    fi

    rank_graft_candidates.py \\
        --asr-dir asr_in \\
        --divergence-tests ${divergence_tests} \\
        \$SEL_ARG \$CONV_ARG \\
        --out-tsv graft_candidate_ranking.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
