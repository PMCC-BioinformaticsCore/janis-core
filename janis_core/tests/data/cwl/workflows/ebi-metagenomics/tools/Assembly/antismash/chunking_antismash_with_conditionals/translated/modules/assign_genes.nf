nextflow.enable.dsl=2

process ASSIGN_GENES {
    debug true
    container "biowardrobe2/iaintersect:v0.0.2"
    publishDir "${params.outdir}/assign_genes"

    input:
    path annotation_filename, stageAs: 'annotation_filename'
    path input_filename, stageAs: 'input_filename'

    output:
    path "*_iaintersect.tsv", emit: result_file

    script:
    """
    iaintersect \
    --a=${annotation_filename} \
    --in=${input_filename} \
    --promoter=${params.promoter_bp} \
    --out=<js>${  let root = inputs.input_filename.basename.split('.').slice(0,-1).join('.');  return (root == "")?inputs.input_filename.basename+"_iaintersect.tsv":root+"_iaintersect.tsv";}</js> \
    """

}
