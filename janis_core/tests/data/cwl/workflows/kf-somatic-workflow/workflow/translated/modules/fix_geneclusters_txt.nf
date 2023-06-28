nextflow.enable.dsl=2

process FIX_GENECLUSTERS_TXT {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/run_antismash/fix_geneclusters_txt"

    input:
    path input_geneclusters_txt, stageAs: 'input_geneclusters_txt'
    val output_filename

    output:
    path "inputs.output_filename", emit: fixed_txt

    script:
    """
    change_geneclusters_ctg.py \
    -i ${input_geneclusters_txt} \
    -o ${output_filename} \
    """

}
