nextflow.enable.dsl=2

process INDEX_STRELKA_BED {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9"
    publishDir "${params.outdir}/index_strelka_bed"

    input:
    path input_file, stageAs: 'input_file'
    path input_index, stageAs: 'input_index'

    output:
    tuple path("*.gz"), path("*.tbi"), emit: output

    script:
    def input_file = input_file ? input_file : ""
    """
    \
    <js>inputs.input_file ? inputs.input_index ? 'echo tabix -p vcf' : 'tabix -p vcf' : 'echo tabix -p vcf'</js> \
    ${input_file} \
    """

}
