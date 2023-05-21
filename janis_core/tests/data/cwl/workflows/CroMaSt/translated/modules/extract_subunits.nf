nextflow.enable.dsl=2

process EXTRACT_SUBUNITS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_subunits"
    cpus "${params.after_qc.rna_prediction.extract_subunits.cpus}"
    memory "${params.after_qc.rna_prediction.extract_subunits.memory}"

    input:
    path input, stageAs: 'input'
    val pattern_lsu
    val pattern_ssu
    val pattern_5.8s
    val pattern_5_s

    output:
    path "sequence-categorisation/*LSU.fasta*", emit: LSU_seqs
    path "sequence-categorisation/*SSU.fasta*", emit: SSU_seqs
    path "sequence-categorisation/*.fa", emit: fastas
    path "sequence-categorisation", optional: true, emit: sequence_categorisation
    stdout, emit: stdout

    script:
    def pattern_5.8s = pattern_5.8s ? "-e ${pattern_5.8s}" : ""
    def pattern_5_s = pattern_5_s ? "-f ${pattern_5_s}" : ""
    """
    get_subunits.py \
    -i ${input} \
    ${pattern_5.8s} \
    ${pattern_5_s} \
    -l ${pattern_lsu} \
    -s ${pattern_ssu} \
    > stdout.txt \
    """

}


process EXTRACT_SUBUNITS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_subunits"
    cpus "${params.after_qc.rna_prediction.extract_subunits.cpus}"
    memory "${params.after_qc.rna_prediction.extract_subunits.memory}"

    input:
    path input, stageAs: 'input'
    val pattern_lsu
    val pattern_ssu
    val pattern_5.8s
    val pattern_5_s

    output:
    path "sequence-categorisation/*LSU.fasta*", emit: LSU_seqs
    path "sequence-categorisation/*SSU.fasta*", emit: SSU_seqs
    path "sequence-categorisation/*.fa", emit: fastas
    path "sequence-categorisation", optional: true, emit: sequence_categorisation
    stdout, emit: stdout

    script:
    def pattern_5.8s = pattern_5.8s ? "-e ${pattern_5.8s}" : ""
    def pattern_5_s = pattern_5_s ? "-f ${pattern_5_s}" : ""
    """
    get_subunits.py \
    -i ${input} \
    ${pattern_5.8s} \
    ${pattern_5_s} \
    -l ${pattern_lsu} \
    -s ${pattern_ssu} \
    > stdout.txt \
    """

}


process EXTRACT_SUBUNITS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/extract_subunits"
    cpus "${params.after_qc.rna_prediction.extract_subunits.cpus}"
    memory "${params.after_qc.rna_prediction.extract_subunits.memory}"

    input:
    path input, stageAs: 'input'
    val pattern_lsu
    val pattern_ssu
    val pattern_5.8s
    val pattern_5_s

    output:
    path "sequence-categorisation/*LSU.fasta*", emit: LSU_seqs
    path "sequence-categorisation/*SSU.fasta*", emit: SSU_seqs
    path "sequence-categorisation/*.fa", emit: fastas
    path "sequence-categorisation", optional: true, emit: sequence_categorisation
    stdout, emit: stdout

    script:
    def pattern_5.8s = pattern_5.8s ? "-e ${pattern_5.8s}" : ""
    def pattern_5_s = pattern_5_s ? "-f ${pattern_5_s}" : ""
    """
    get_subunits.py \
    -i ${input} \
    ${pattern_5.8s} \
    ${pattern_5_s} \
    -l ${pattern_lsu} \
    -s ${pattern_ssu} \
    > stdout.txt \
    """

}
