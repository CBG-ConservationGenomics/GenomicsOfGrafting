process ANNOTATE_DEEPTMHMM {
    label 'process_medium'
    publishDir "${params.outdir}/04_annotation/deeptmhmm", mode: 'copy'

    input:
    path fasta

    output:
    path "deeptmhmm.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    params.run_deeptmhmm

    script:
    """
    set -euo pipefail

    # DeepTMHMM is commonly run via BioLib. Override params.deeptmhmm_bin on your server.
    # Expected: produces predicted_topologies.3line or similar.
    mkdir -p deeptmhmm_out

    if command -v biolib >/dev/null 2>&1; then
        biolib run DTU/DeepTMHMM --fasta ${fasta} || true
        # BioLib often writes into a biolib results directory
        pred=\$(find . -name 'predicted_topologies.3line' -o -name '*deeptmhmm*' -type f | head -n1 || true)
        if [ -n "\$pred" ]; then
            cp "\$pred" deeptmhmm_out/predicted_topologies.3line
        fi
    fi

    if [ -f deeptmhmm_out/predicted_topologies.3line ]; then
        # Convert 3line to TSV: id, topology_line
        python3 - <<'PY'
from pathlib import Path
lines = Path("deeptmhmm_out/predicted_topologies.3line").read_text().splitlines()
out = open("deeptmhmm.tsv", "w")
out.write("id\\ttopology\\n")
i = 0
while i < len(lines):
    if lines[i].startswith(">"):
        pid = lines[i][1:].split()[0]
        topo = lines[i+2] if i+2 < len(lines) else ""
        out.write(f"{pid}\\t{topo}\\n")
        i += 3
    else:
        i += 1
out.close()
PY
    else
        echo -e "id\\ttopology" > deeptmhmm.tsv
        echo "WARNING: DeepTMHMM output not found; wrote empty table. Configure biolib/DeepTMHMM on the server." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptmhmm: biolib/DeepTMHMM
    END_VERSIONS
    """
}
