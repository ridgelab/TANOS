process JACKKNIFE_SUPPORT_COUNT {

    cpus 2
    memory '8 GB'
    time '4d'

    input:
    tuple val(replicate_count), path(replicate_treefiles), path(main_tree), path(support_script)

    output:
    tuple val(replicate_count), path("${replicate_count}_replicates.treefile"), emit: support_tree

    script:
    """
    set -euo pipefail

    cat ${replicate_treefiles} > jackknife_replicates.trees
    python3 "${support_script}" \\
        "${main_tree}" \\
        jackknife_replicates.trees \\
        "${replicate_count}_replicates.treefile"
    rm jackknife_replicates.trees
    """
}
