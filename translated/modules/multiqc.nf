nextflow.enable.dsl=2

process MULTIQC {
    
    container "quay.io/biocontainers/multiqc:1.7--py_4"
    publishDir "${params.outdir}/multiqc"

    input:
    path unknown1
    path unknown2

    output:
    path "report.html", emit: outHtmlReport
    path "unknown", emit: outStats

    script:
    """
    multiqc multiqc_WDir \
    """

}
