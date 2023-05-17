nextflow.enable.dsl=2

process FOLDER_COMPRESSION {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/folder_compression"

    input:
    path indir, stageAs: 'indir'

    output:
    path "{inputs.indir.basename}.tar.gz", emit: outfile

    script:
    """
    tar czfh \
    """

}
