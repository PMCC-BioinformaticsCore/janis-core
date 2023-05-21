nextflow.enable.dsl=2

process DIAMOND_RUN {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/diamond/diamond_run"
    cpus "${params.after_qc.functional_annotation_and_post_processing.diamond.diamond_run.cpus}"
    memory "${params.after_qc.functional_annotation_and_post_processing.diamond.diamond_run.memory}"

    input:
    path query_input_file, stageAs: 'query_input_file'
    val database_file
    val max_target_seqs
    val output_format
    val strand
    val threads

    output:
    path "{inputs.queryInputFile.basename}.diamond_matches", emit: matches

    script:
    def max_target_seqs = max_target_seqs ? "--max-target-seqs ${max_target_seqs}" : ""
    def output_format = output_format ? "--outfmt ${output_format}" : ""
    def strand = strand ? "--strand ${strand}" : ""
    def threads = threads ? "--threads ${threads}" : ""
    """
    diamond blastp \
    ${strand} \
    --query ${query_input_file} \
    --db ${database_file} \
    ${max_target_seqs} \
    ${output_format} \
    ${threads} \
    --out "${query_input_file.name}.diamond_matches" \
    """

}
