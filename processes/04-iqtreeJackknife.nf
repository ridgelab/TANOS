process IQTREE_JACKKNIFE {

    tag "${taxon}-${rep}"

    cpus 16
    memory '24 GB'
    time '1d'

    publishDir "${params.outdir}/${taxon}", mode: 'copy'

    input:
    tuple path(fa), val(rep)

    output:
    path "${taxon}-${rep}.*"

    script:
    """
    iqtree2 \
        -nt ${task.cpus} \
        -s "${aln}" \
        -t RANDOM \
        -pre "${taxon}-${rep}" \
        -m "${params.model}"
    """
}