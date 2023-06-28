nextflow.enable.dsl=2

process BAM_INDEX {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/bam_index"

    input:
    path bam_file, stageAs: 'bam_file'

    output:
    path "{inputs.bam_file.basename}.bai", emit: bam_index

    script:
    def threads = params.threads ? params.threads : 2
    """
    samtools index \
    -@ ${threads} \
    ${bam_file} \
    "${bam_file.name}.bai" \
    """

}
