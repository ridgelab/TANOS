process JACKKNIFE_SUPPORT {

    tag "${replicate_count} replicates"

    cpus 4
    memory '8 GB'
    time '1d'

    publishDir "${params.outdir}/jackknifeSupport", mode: 'copy'

    input:
    tuple val(replicate_count), path(jackknife_treefiles)
    path main_tree

    output:
    tuple val(replicate_count), path("${replicate_count}_replicates.suptree"), emit: support_tree

    script:
    """
    cat ${jackknife_treefiles} > jackknife_replicates.trees

    iqtree2 \\
        -nt ${task.cpus} \\
        -pre ${replicate_count}_replicates \\
        -z jackknife_replicates.trees \\
        -t ${main_tree} \\
        -sup ${main_tree}
    """
}
