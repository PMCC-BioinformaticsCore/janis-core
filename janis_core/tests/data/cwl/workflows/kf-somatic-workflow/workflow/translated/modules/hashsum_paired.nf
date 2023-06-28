nextflow.enable.dsl=2

process HASHSUM_PAIRED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/hashsum_paired"

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


process HASHSUM_PAIRED {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/hashsum_paired"

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
