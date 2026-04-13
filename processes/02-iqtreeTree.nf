process MAIN_TREE {

    cpus 16
    memory '24 GB'
    time '1d'

    publishDir "${params.outdir}/mainTree", mode: 'copy'

    input:
    path input_aln
    val model

    output:
    path "tree.*", emit: files

    script:
    """
    set -e

    # threads + memory from Nextflow
    THREADS=${task.cpus}
    MEM_GB=${task.memory.toGiga()}
    OUTPUT_PFX='mainTree'

    # run iqtree
    iqtree2 \
        -nt \${THREADS} \
        -mem \${MEM_GB}G \
        -s "${input_aln}" \
        -pre \${OUTPUT_PFX} \
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