nextflow.enable.dsl=2

process SAMTOOLS_FAIDX {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"
    publishDir "${params.outdir}/prepare_reference/samtools_faidx"

    input:
    path input_fasta, stageAs: 'input_fasta'
    path input_index, stageAs: 'input_index'

    output:
    path "*.fai", emit: fai

    script:
    """
    \
    <js>inputs.input_index ? 'echo samtools faidx' : 'samtools faidx' </js> \
    ${input_fasta} \
    """

}
