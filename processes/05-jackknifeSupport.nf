process JACKKNIFE_SUPPORT {

    cpus 4
    memory '8 GB'
    time '1d'

    publishDir "${params.outdir}/jackknifeSupport", mode: 'copy'

    input:
    path input_aln
    path jackknife_treefiles
    path main_tree
    val model

    output:
    path "primary_with_jackknife.*"

    script:
    """
        cat ${jackknife_treefiles} > jackknife_replicates.trees
        
        iqtree2 \
            -nt ${task.cpus} \
            -s ${input_aln} \
            -pre primary_with_jackknife \
            -z jackknife_replicates.trees \
            -t ${main_tree} \
            -sup ${main_tree} \
            -m "${model}"
    """
}
