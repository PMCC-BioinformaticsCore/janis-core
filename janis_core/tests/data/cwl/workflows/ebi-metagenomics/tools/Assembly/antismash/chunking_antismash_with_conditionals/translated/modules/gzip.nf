nextflow.enable.dsl=2

process GZIP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/gzip"

    input:
    path infile, stageAs: 'infile'

    output:
    path "{inputs.infile.basename}.gz", emit: output

    script:
    """
    gzip -c \
    > "${infile.name}.gz" \
    """

}
