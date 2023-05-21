nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process MOVE_TO_SEQ_CAT_FOLDER {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/move_to_seq_cat_folder"

    input:
    path file_list

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['out']}", emit: out

    script:
    """
    nodejs return_directory.js \
    > cwl.output.json \
    """

}


import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process MOVE_TO_SEQ_CAT_FOLDER {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/move_to_seq_cat_folder"

    input:
    path file_list

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['out']}", emit: out

    script:
    """
    nodejs return_directory.js \
    > cwl.output.json \
    """

}
