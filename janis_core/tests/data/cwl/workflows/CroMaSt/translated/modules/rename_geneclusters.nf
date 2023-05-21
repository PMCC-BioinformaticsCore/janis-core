nextflow.enable.dsl=2

process RENAME_GENECLUSTERS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/rename_geneclusters"
    memory "${params.after_qc.antismash.chunking.rename_geneclusters.memory}"

    input:
    path initial_file, stageAs: 'initial_file'
    val out_file_name

    output:
    path "inputs.out_file_name", emit: renamed_file

    script:
    """
    mv \
    ${initial_file} \
    ${out_file_name} \
    """

}
