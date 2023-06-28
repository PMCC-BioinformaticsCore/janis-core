nextflow.enable.dsl=2

process COMPRESS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/compress"

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
