nextflow.enable.dsl=2

process FIX_EMBL_AND_GBK {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/run_antismash/fix_embl_and_gbk"
    memory "${params.after_qc.antismash.chunking.run_antismash.fix_embl_and_gbk.memory}"

    input:
    path embl_file, stageAs: 'embl_file'
    path names_table, stageAs: 'names_table'
    val embl_filename
    val gbk_filename

    output:
    path "inputs.embl_filename", emit: fixed_embl
    path "inputs.gbk_filename", emit: fixed_gbk

    script:
    """
    change_antismash_output.py \
    -i ${embl_file} \
    -t ${names_table} \
    -e ${embl_filename} \
    -g ${gbk_filename} \
    """

}
