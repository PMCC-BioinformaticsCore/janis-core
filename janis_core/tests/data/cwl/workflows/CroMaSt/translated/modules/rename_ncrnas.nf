nextflow.enable.dsl=2

process RENAME_NCRNAS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/rename_ncrnas"
    memory "${params.after_qc.other_ncrnas.rename_ncrnas.memory}"

    input:
    path initial_file, stageAs: 'initial_file'

    output:
    path "inputs.out_file_name", emit: renamed_file

    script:
    """
    mv \
    ${initial_file} \
    ${params.after_qc.other_ncrnas.rename_ncrnas_out_file_name} \
    """

}


process RENAME_NCRNAS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/rename_ncrnas"
    memory "${params.after_qc.other_ncrnas.rename_ncrnas.memory}"

    input:
    path initial_file, stageAs: 'initial_file'

    output:
    path "inputs.out_file_name", emit: renamed_file

    script:
    """
    mv \
    ${initial_file} \
    ${params.after_qc.other_ncrnas.rename_ncrnas_out_file_name} \
    """

}
