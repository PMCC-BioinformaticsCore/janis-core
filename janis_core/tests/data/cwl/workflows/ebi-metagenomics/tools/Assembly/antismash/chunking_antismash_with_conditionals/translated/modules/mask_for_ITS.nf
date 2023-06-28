nextflow.enable.dsl=2

process MASK_FOR_ITS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/its/mask_for_its"
    cpus "${params.after_qc.its.mask_for_its.cpus}"
    memory "${params.after_qc.its.mask_for_its.memory}"

    input:
    path maskfile, stageAs: 'maskfile'
    path sequences, stageAs: 'sequences'

    output:
    path "ITS_masked.fasta", emit: masked_sequences

    script:
    """
    bedtools maskfasta \
    -fo ITS_masked.fasta \
    -fi ${sequences} \
    -bed ${maskfile} \
    """

}
