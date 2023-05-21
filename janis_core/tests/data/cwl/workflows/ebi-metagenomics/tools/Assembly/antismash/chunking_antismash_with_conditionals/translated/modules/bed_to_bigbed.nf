nextflow.enable.dsl=2

process BED_TO_BIGBED {
    debug true
    container "biowardrobe2/ucscuserapps:v358"
    publishDir "${params.outdir}/bed_to_bigbed"

    input:
    path chrom_length_file, stageAs: 'chrom_length_file'
    path input_bed, stageAs: 'input_bed'
    val output_filename

    output:
    path "*", emit: bigbed_file

    script:
    """
    bedToBigBed \
    -type=${params.bed_to_bigbed_bed_type} \
    ${input_bed} \
    ${chrom_length_file} \
    ${output_filename} \
    """

}
