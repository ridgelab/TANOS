process MODEL_TEST {

    tag "${job_name}"

    cpus 16
    memory '24 GB'
    time '1d'

    publishDir "results/modelTest", mode: 'copy'

    input:
    val job_name
    path input_aln

    output:
    path "modelTest.*", emit: files
    path "model.txt", emit: model

    script:
    """
    set -e

    # threads + memory from Nextflow
    THREADS=${task.cpus}
    MEM_GB=${task.memory.toGiga()}
    OUTPUT_PFX='results/modelTest'

    # run iqtree
    iqtree2 \
        -nt \${THREADS} \
        -mem \${MEM_GB}G \
        -s "${input_aln}" \
        -pre \${OUTPUT_PFX} \
        -m TESTONLY

    # handle outputs like your script
    if [ \$? -eq 0 ]; then
        chmod 444 \${OUTPUT_PFX}.*
    else
        rm -f \${OUTPUT_PFX}.*
        exit 1
    fi

    # Extract best model
    awk -F': ' '/Best-fit model according to BIC:/ {print \$2}' \
        \${OUTPUT_PFX}.iqtree > model.txt
    """
}