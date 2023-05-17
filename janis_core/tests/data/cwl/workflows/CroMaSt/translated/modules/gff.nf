nextflow.enable.dsl=2

process GFF {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/gff"
    memory "${params.after_qc.functional_annotation_and_post_processing.gff.memory}"

    input:
    path eggnog_results, stageAs: 'eggnog_results'
    path input_faa, stageAs: 'input_faa'
    path ips_results, stageAs: 'ips_results'
    val output_name

    output:
    path "{inputs.output_name}.bgz", emit: output_gff_gz
    path "{inputs.output_name}.bgz.tbi", emit: output_gff_index
    stdout, emit: stdout

    script:
    """
    build_assembly_gff.py \
    -e ${eggnog_results} \
    -f ${input_faa} \
    -i ${ips_results} \
    -o ${output_name} \
    > stdout.txt \
    """

}
