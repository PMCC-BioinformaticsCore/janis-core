nextflow.enable.dsl=2

process ANTISMASH_GFF {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/antismash/chunking/antismash_gff"

    input:
    path antismash_embl, stageAs: 'antismash_embl'
    path antismash_geneclus, stageAs: 'antismash_geneclus'
    val output_name

    output:
    path "{inputs.output_name}.bgz", emit: output_gff_bgz
    path "{inputs.output_name}.bgz.tbi", emit: output_gff_index
    stdout, emit: stdout

    script:
    """
    antismash_to_gff.py \
    -e ${antismash_embl} \
    -g ${antismash_geneclus} \
    -o ${output_name} \
    > stdout.txt \
    """

}
