nextflow.enable.dsl=2

process SUPPRESS_TAX {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/suppress_tax"
    memory "${params.after_qc.suppress_tax.memory}"

    input:
    path its_file, stageAs: 'its_file'
    path lsu_file, stageAs: 'lsu_file'
    path ssu_file, stageAs: 'ssu_file'
    path its_dir, stageAs: 'its_dir'
    path lsu_dir, stageAs: 'lsu_dir'
    path ssu_dir, stageAs: 'ssu_dir'

    output:
    stdout, emit: its_length
    path "*.fasta.gz", optional: true, emit: out_fastas_tax
    path "suppressed", optional: true, emit: out_suppress
    path "taxonomy-summary", optional: true, emit: out_tax

    script:
    def its_dir = its_dir ? "--its-dir ${its_dir}" : ""
    def its_file = its_file ? "--its-file ${its_file}" : ""
    def lsu_dir = lsu_dir ? "--lsu-dir ${lsu_dir}" : ""
    def lsu_file = lsu_file ? "--lsu-file ${lsu_file}" : ""
    def ssu_dir = ssu_dir ? "--ssu-dir ${ssu_dir}" : ""
    def ssu_file = ssu_file ? "--ssu-file ${ssu_file}" : ""
    """
    its-length-new.py \
    ${its_file} \
    ${lsu_file} \
    ${ssu_file} \
    ${its_dir} \
    ${lsu_dir} \
    ${ssu_dir} \
    > ITS_LENGTH \
    """

}
