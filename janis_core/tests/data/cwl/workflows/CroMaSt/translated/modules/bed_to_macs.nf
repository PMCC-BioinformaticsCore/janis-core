nextflow.enable.dsl=2

process BED_TO_MACS {
    debug true
    container "biowardrobe2/scidap:v0.0.3"
    publishDir "${params.outdir}/bed_to_macs"

    input:
    path input_file, stageAs: 'input_file'

    output:
    path "*", emit: output_file

    script:
    """
    bash -c \
    cat $0 | grep -v "#" | awk 'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"} {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t0\t0\t0\t0\t"$4}' > `basename $0` \
    ${input_file} \
    """

}
