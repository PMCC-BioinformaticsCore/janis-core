nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process RETURN_ITS_DIR {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/return_its_dir"

    input:
    path dir_list

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['out']}", emit: out

    script:
    """
    nodejs return_directory.js \
    > cwl.output.json \
    """

}
