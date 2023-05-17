nextflow.enable.dsl=2

process MAKE_TAB_SEP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/make_tab_sep"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.run_hmmer.make_tab_sep.memory}"

    input:
    path input_table, stageAs: 'input_table'

    output:
    path "*.tsv", emit: output_with_tabs

    script:
    """
    hmmscan_tab.py \
    -i ${input_table} \
    -o "${${input_table}.baseName}.tsv" \
    """

}


process MAKE_TAB_SEP {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/functional_annotation_and_post_processing/functional_annotation/run_hmmer/make_tab_sep"
    memory "${params.after_qc.functional_annotation_and_post_processing.functional_annotation.run_hmmer.make_tab_sep.memory}"

    input:
    path input_table, stageAs: 'input_table'

    output:
    path "*.tsv", emit: output_with_tabs

    script:
    """
    hmmscan_tab.py \
    -i ${input_table} \
    -o "${${input_table}.baseName}.tsv" \
    """

}
