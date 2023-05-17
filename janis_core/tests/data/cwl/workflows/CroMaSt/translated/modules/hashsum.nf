nextflow.enable.dsl=2

process HASHSUM {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/hashsum"
    cpus "${params.before_qc.hashsum.cpus}"
    memory "${params.before_qc.hashsum.memory}"

    input:
    path input_file, stageAs: 'input_file'

    output:
    path "*sha1", emit: hashsum

    script:
    def input_file = input_file ? "-i ${input_file}" : ""
    """
    generate_checksum.py \
    ${input_file} \
    """

}
