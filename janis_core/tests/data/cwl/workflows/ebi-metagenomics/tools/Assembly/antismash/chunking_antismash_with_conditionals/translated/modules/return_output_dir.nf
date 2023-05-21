nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process RETURN_OUTPUT_DIR {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/return_output_dir"

    input:
    path file_list
    val dir_name

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

process RETURN_OUTPUT_DIR {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/return_output_dir"

    input:
    path file_list
    val dir_name

    output:
    path "${jsonSlurper.parseText(file("${task.workDir}/cwl.output.json").text)['out']}", emit: out

    script:
    """
    nodejs return_directory.js \
    > cwl.output.json \
    """

}
