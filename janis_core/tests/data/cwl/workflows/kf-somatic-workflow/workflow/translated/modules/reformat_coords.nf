nextflow.enable.dsl=2

process REFORMAT_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/reformat_coords"

    input:
    path all_coordinates, stageAs: 'all_coordinates'

    output:
    stdout, emit: maskfile

    script:
    """
    format_bedfile \
    -i ${all_coordinates} \
    > ITS-maskfile \
    """

}
