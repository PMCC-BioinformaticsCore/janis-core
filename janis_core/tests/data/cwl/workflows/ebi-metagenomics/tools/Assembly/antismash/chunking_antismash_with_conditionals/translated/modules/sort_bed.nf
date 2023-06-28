nextflow.enable.dsl=2

process SORT_BED {
    debug true
    container "biowardrobe2/scidap:v0.0.3"
    publishDir "${params.outdir}/sort_bed"

    input:
    path unsorted_file, stageAs: 'unsorted_file'

    output:
    path "*", emit: sorted_file

    script:
    def key = params.sort_bed_key.join(' ')
    """
    sort \
    ${key} \
    ${unsorted_file} \
    > stdout.txt \
    """

}
