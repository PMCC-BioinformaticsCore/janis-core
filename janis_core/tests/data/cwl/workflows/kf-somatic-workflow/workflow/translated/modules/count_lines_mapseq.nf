nextflow.enable.dsl=2

process COUNT_LINES_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/count_lines_mapseq"

    input:
    path input_file, stageAs: 'input_file'

    output:
    val "count", emit: number

    script:
    """
    bash \
    -c "expr \\<js>cat $(inputs.input_file.path) | wc -l</js>
    " \
    > count \
    """

}


process COUNT_LINES_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/count_lines_mapseq"

    input:
    path input_file, stageAs: 'input_file'

    output:
    val "count", emit: number

    script:
    """
    bash \
    -c "expr \\<js>cat $(inputs.input_file.path) | wc -l</js>
    " \
    > count \
    """

}


process COUNT_LINES_MAPSEQ {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/count_lines_mapseq"

    input:
    path input_file, stageAs: 'input_file'

    output:
    val "count", emit: number

    script:
    """
    bash \
    -c "expr \\<js>cat $(inputs.input_file.path) | wc -l</js>
    " \
    > count \
    """

}
