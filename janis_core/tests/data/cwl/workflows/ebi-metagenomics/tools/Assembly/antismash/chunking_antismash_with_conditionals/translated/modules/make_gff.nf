nextflow.enable.dsl=2

process MAKE_GFF {
    debug true
    container "biowardrobe2/rose:v0.0.2"
    publishDir "${params.outdir}/make_gff"

    input:
    path islands_file, stageAs: 'islands_file'
    path islands_control_file, stageAs: 'islands_control_file'

    output:
    path "*", emit: gff_file

    script:
    def islands_control_file = islands_control_file ? islands_control_file : ""
    """
    makegff \
    ${islands_file} \
    <js>${ let root = inputs.islands_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.islands_file.basename+".gff":root+".gff"; }</js> \
    ${islands_control_file} \
    """

}
