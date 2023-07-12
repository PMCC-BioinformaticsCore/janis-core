nextflow.enable.dsl=2

process ADD_ISLAND_NAMES {
    debug true
    container "biowardrobe2/scidap:v0.0.3"
    publishDir "${params.outdir}/add_island_names"

    input:
    path input_file
    val param

    output:
    path "*", emit: output_file

    script:
    def input_file = input_file.join(' ')
    """
    bash -c \
    echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tchrom\tstart\tend\tlength\tregion\tname\tscore" > `basename $2`; cat $0 | grep -v refseq_id | paste - $1 | cut -f 1-9,15,19,20 >> `basename $2` \
    ${input_file} \
    ${param} \
    """

}
