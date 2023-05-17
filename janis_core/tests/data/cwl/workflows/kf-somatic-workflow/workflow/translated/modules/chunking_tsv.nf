nextflow.enable.dsl=2

process CHUNKING_TSV {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/chunking_tsv"

    input:
    path infile

    output:
    path "{inputs.outdirname}/*", emit: chunks

    script:
    def infile = infile.join(' ')
    """
    run_result_file_chunker.py \
    -i ${infile} \
    -f ${params.after_qc.functional_annotation_and_post_processing.chunking_tsv_format_file} \
    -o ${params.after_qc.functional_annotation_and_post_processing.chunking_tsv_outdirname} \
    """

}
