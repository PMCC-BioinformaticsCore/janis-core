nextflow.enable.dsl=2

process VCF_COMPRESS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/vcf_compress"

    input:
    path inputfile, stageAs: 'inputfile'

    output:
    path "{inputs.inputfile.basename}.gz", emit: outfile

    script:
    def threads = params.threads ? params.threads : 1
    """
    pigz -c \
    -p ${threads} \
    > "${inputfile.name}.gz" \
    """

}
