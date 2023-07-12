nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process EXPRESSIONTOOL_BAM_INDEX {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/expressiontool_bam_index"

    input:
    path bam_file, stageAs: 'bam_file'
    path bam_index, stageAs: 'bam_index'

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['hybrid_bamindex']}", emit: hybrid_bamindex

    script:
    """
    nodejs expression_bam_index.js \
    > cwl.output.json \
    """

}
