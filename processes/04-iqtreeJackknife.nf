process IQTREE_JACKKNIFE {

    tag "${fa.simpleName}-${rep}"

    cpus 4
    memory '8 GB'
    time '1d'

    publishDir "${params.outdir}/${fa.simpleName}", mode: 'copy'

    input:
    tuple path(fa), val(rep)
    val model

    output:
    path "${fa.simpleName}-${rep}.*"

    script:
    """
    iqtree2 \
        -nt ${task.cpus} \
        -mem ${task.memory.toGiga()}G \
        -s "${fa}" \
        -t RANDOM \
        -pre "${fa.simpleName}-${rep}" \
        -m "${model}"
    """
}