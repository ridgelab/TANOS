process JACKKNIFE_CONVERGENCE {

    cpus 1
    memory '2 GB'
    time '1h'

    publishDir "${params.outdir}/jackknifeSupport", mode: 'copy'

    input:
    path support_trees
    path analyzer_script

    output:
    stdout emit: summary
    path "jackknife_convergence.csv"
    path "jackknife_convergence.png"

    script:
    """
    set -euo pipefail
    mkdir -p convergence_results

    for support_tree in *_replicates.treefile; do
        replicate_count="\${support_tree%%_*}"
        result_dir="convergence_results/\${replicate_count}_replicates"
        mkdir -p "\$result_dir"
        ln -s "../../\$support_tree" \\
            "\$result_dir/primary_with_support.treefile"
    done

    python3 "${analyzer_script}" convergence_results \\
        --csv jackknife_convergence.csv \\
        --png jackknife_convergence.png
    rm -rf convergence_results
    """
}
