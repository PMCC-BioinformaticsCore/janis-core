nextflow.enable.dsl=2

process NANOPLOT {
    
    container "quay.io/biocontainers/nanoplot:1.28.2--py_0"
    publishDir "${params.outdir}/nanoplot"

    input:
    path unknown1

    output:
    path "LogTransformed_HistogramReadlength.*", emit: outLogReadLength
    path "NanoStats.txt", emit: outNanostats
    path "NanoStats_post_filtering.txt", emit: outNanostatsPostFiltering
    path "HistogramReadlength.*", emit: outReadLength
    path "NanoPlot-report.html", emit: outputHtml

    script:
    """
    NanoPlot \
    """

}
