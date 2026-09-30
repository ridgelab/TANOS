process IQTREE_JACKKNIFE {

    tag "${fa.simpleName}-${rep}"

    cpus 8
    memory '16 GB'
    time '4d'

    publishDir "${params.outdir}/${fa.simpleName}", mode: 'copy'

    input:
    tuple path(fa), val(rep)
    val model

    output:
    path "${fa.simpleName}-${rep}.treefile"

    script:
    """
    iqtree2 \
        -nt ${task.cpus} \
        -mem ${task.memory.toGiga()}G \
        -s "${fa}" \
        -t RANDOM \
        -pre "${fa.simpleName}-${rep}" \
        -m "${model}" \
        -fast
    """
}
