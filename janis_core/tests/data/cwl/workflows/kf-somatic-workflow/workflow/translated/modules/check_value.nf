nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process CHECK_VALUE {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/check_value"

    input:
    val number

    output:
    val "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['out']}", emit: out

    script:
    """
    nodejs check_value.js \
    > cwl.output.json \
    """

}
