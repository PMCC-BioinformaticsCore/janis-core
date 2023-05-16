nextflow.enable.dsl=2

process QUAST {
    
    container "ppp-janis-translate:quast-5.0.2"
    publishDir "${params.outdir}/quast"

    input:
    path unknown1

    output:
    path "outputdir/quast.log", emit: outLog
    path "outputdir/report.tsv", emit: outQuastTabular
    path "outputdir/report.html", emit: outReportHtml
    path "outputdir/report.pdf", emit: outReportPdf

    script:
    """
    quast metaquast \
    """

}
