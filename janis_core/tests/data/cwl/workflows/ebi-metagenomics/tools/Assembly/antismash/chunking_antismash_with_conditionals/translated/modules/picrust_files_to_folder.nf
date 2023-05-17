nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process PICRUST_FILES_TO_FOLDER {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/picrust_files_to_folder"

    input:
    path files
    path folders

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['results']}", emit: results

    script:
    """
    nodejs files_to_folder.js \
    > cwl.output.json \
    """

}
