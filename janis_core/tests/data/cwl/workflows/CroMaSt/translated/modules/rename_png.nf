nextflow.enable.dsl=2

process RENAME_PNG {
    debug true
    container "biowardrobe2/scidap:v0.0.3"
    publishDir "${params.outdir}/rename_png"

    input:
    path source_file, stageAs: 'source_file'
    val target_filename

    output:
    path "*", emit: target_file

    script:
    """
    cp \
    ${source_file} \
    ${target_filename} \
    """

}
