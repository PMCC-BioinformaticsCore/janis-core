nextflow.enable.dsl=2

process CARVEME {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/carveme"

    input:
    path protein_file, stageAs: 'protein_file'
    path mediadb, stageAs: 'mediadb'

    output:
    path "{inputs.protein_file.nameroot}.GEM.xml", optional: true, emit: carveme_gem

    script:
    def gapfill = params.gapfill ? "--gapfill ${params.gapfill}" : ""
    def mediadb = mediadb ? "--mediadb ${mediadb}" : ""
    """
    carve \
    ${mediadb} \
    ${gapfill} \
    --output "${${protein_file}.baseName}.GEM.xml" \
    ${protein_file} \
    --fbc2 \
    """

}
