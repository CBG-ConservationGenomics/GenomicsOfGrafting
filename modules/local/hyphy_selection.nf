process HYPHY_SELECTION {
    tag "${og_id}"
    label 'process_high'
    publishDir "${params.outdir}/10_graft_phenotype/selection/${og_id}", mode: 'copy'

    input:
    tuple val(og_id), path(codon_alignment), path(tree)
    path foreground_tips

    output:
    tuple val(og_id), path("${og_id}.hyphy_summary.tsv"), emit: summary
    path "*.json", emit: json, optional: true
    path "versions.yml", emit: versions

    when:
    params.run_hyphy

    script:
    """
    set -euo pipefail

    # Codon alignment + tree required. Foreground tips labeled for aBSREL/RELAX via trees.
    # Build a simple labeled tree copy: tips in foreground_tips get {Foreground}.
    python3 - <<'PY'
from Bio import Phylo
from pathlib import Path
fg = {ln.strip() for ln in open("${foreground_tips}") if ln.strip()}
tree = Phylo.read("${tree}", "newick")
for tip in tree.get_terminals():
    name = tip.name.split("|")[0] if tip.name and "|" in tip.name else tip.name
    if name in fg:
        tip.name = f"{tip.name}{{Foreground}}"
Phylo.write(tree, "${og_id}.fg.nwk", "newick")
PY

    echo -e "Orthogroup\\tmethod\\tpvalue\\tkey_stat\\tstatus" > ${og_id}.hyphy_summary.tsv

    run_one() {
        local method=\$1
        local outjson=\$2
        if hyphy \$method --alignment ${codon_alignment} --tree ${og_id}.fg.nwk CPU=${task.cpus} > ${method}.log 2>&1; then
            # Best-effort scrape; HyPhy JSON names vary by method
            pv=\$(python3 - <<PY
import json,glob,re
files=glob.glob("*.json")+glob.glob("${outjson}")
pval=""; stat=""
for f in files:
    try:
        d=json.load(open(f))
    except Exception:
        continue
    txt=json.dumps(d)
    m=re.search(r'"p[^"]*[Pp](?:alue)?"\\s*:\\s*([0-9.eE+-]+)', txt)
    if m: pval=m.group(1); break
print(pval)
PY
)
            echo -e "${og_id}\\t\${method}\\t\${pv}\\t\\tok" >> ${og_id}.hyphy_summary.tsv
        else
            echo -e "${og_id}\\t\${method}\\t\\t\\tfailed" >> ${og_id}.hyphy_summary.tsv
        fi
    }

    # Episodic selection / gene-wide / site-level / relaxation
    run_one BUSTED ${og_id}.BUSTED.json || true
    run_one aBSREL ${og_id}.aBSREL.json || true
    run_one MEME ${og_id}.MEME.json || true
    run_one RELAX ${og_id}.RELAX.json || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hyphy: \$(hyphy --version 2>&1 | head -n1 || echo NA)
    END_VERSIONS
    """
}
