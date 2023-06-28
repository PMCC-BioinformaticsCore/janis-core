nextflow.enable.dsl=2

process EXTRACT_SUBUNITS_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_subunits_coords"
    cpus "${params.after_qc.rna_prediction.extract_subunits_coords.cpus}"
    memory "${params.after_qc.rna_prediction.extract_subunits_coords.memory}"

    input:
    path input, stageAs: 'input'
    val pattern_lsu
    val pattern_ssu

    output:
    path "*LSU*", emit: LSU_seqs
    path "*SSU*", emit: SSU_seqs
    path "RNA-counts", emit: counts
    stdout, emit: stdout

    script:
    """
    get_subunits_coords.py \
    -i ${input} \
    -l ${pattern_lsu} \
    -s ${pattern_ssu} \
    > stdout.txt \
    """

}


process EXTRACT_SUBUNITS_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_subunits_coords"
    cpus "${params.after_qc.rna_prediction.extract_subunits_coords.cpus}"
    memory "${params.after_qc.rna_prediction.extract_subunits_coords.memory}"

    input:
    path input, stageAs: 'input'
    val pattern_lsu
    val pattern_ssu

    output:
    path "*LSU*", emit: LSU_seqs
    path "*SSU*", emit: SSU_seqs
    path "RNA-counts", emit: counts
    stdout, emit: stdout

    script:
    """
    get_subunits_coords.py \
    -i ${input} \
    -l ${pattern_lsu} \
    -s ${pattern_ssu} \
    > stdout.txt \
    """

}
