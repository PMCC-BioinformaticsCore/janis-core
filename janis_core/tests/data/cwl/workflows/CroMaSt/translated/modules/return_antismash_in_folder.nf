nextflow.enable.dsl=2

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

process RETURN_ANTISMASH_IN_FOLDER {
    debug true
    container "node:latest"
    publishDir "${params.outdir}/after_qc/antismash/no_antismash_subwf/return_antismash_in_folder"

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
