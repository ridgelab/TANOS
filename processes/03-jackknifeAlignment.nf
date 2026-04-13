process JACKKNIFE_ALIGNMENT {

    cpus 1
    memory '4 GB'
    time '1h'

    publishDir "${params.outdir}/jackknife/tree", mode: 'copy'

    input:
    path input_aln

    output:
    path "*.phy", emit: files

    script:
    """
    set -e

    python3 ${projectDir}/scripts/03-jackknifeAlignment.py \
        "${input_aln}" \
        .
    """
}