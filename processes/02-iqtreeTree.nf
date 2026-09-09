process MAIN_TREE {

    cpus 16
    memory '24 GB'
    time '1d'

    publishDir "${params.outdir}/mainTree", mode: 'copy'

    input:
    path input_aln
    val model

    output:
    path "mainTree.*", emit: files

    script:
    """
    set -e

    # run iqtree
    iqtree2 \
        -nt ${task.cpus} \
        -mem ${task.memory.toGiga()}G \
        -s "${input_aln}" \
        -pre 'mainTree' \
        -m "${model}" 

    # handle outputs like your script
    if [ \$? -eq 0 ]; then
        chmod 444 \${OUTPUT_PFX}.*
    else
        rm -f \${OUTPUT_PFX}.*
        exit 1
    fi
    """
}