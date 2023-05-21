nextflow.enable.dsl=2

process RUN_ANTISMASH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/run_antismash/run_antismash"

    input:
    path input_fasta, stageAs: 'input_fasta'
    val accession
    val chunk_num

    output:
    path "{inputs.outdirname}/*final.embl", emit: embl_file
    path "{inputs.outdirname}/*final.gbk", emit: gbk_file
    path "{inputs.outdirname}/geneclusters.js", emit: geneclusters_js
    path "{inputs.outdirname}/geneclusters.txt", emit: geneclusters_txt

    script:
    """
    run_antismash_short.sh \
    -i ${input_fasta} \
    -o ${params.after_qc.antismash.chunking.run_antismash.run_antismash_outdirname} \
    2> stderr.txt \
    > stdout.txt \
    """

}
