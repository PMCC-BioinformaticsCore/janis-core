nextflow.enable.dsl=2

process REDUCE_BED {
    debug true
    container "biowardrobe2/scidap:v0.0.3"
    publishDir "${params.outdir}/reduce_bed"

    input:
    path input_file, stageAs: 'input_file'

    output:
    path "*", emit: output_file

    script:
    """
    bash -c \
    cat $0 | cut -f 1-4 > `basename $0` \
    ${input_file} \
    """

}
