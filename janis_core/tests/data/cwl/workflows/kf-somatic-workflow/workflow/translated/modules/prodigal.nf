nextflow.enable.dsl=2

process PRODIGAL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/prodigal"

    input:
    path input_fasta, stageAs: 'input_fasta'

    output:
    path "{inputs.input_fasta.nameroot}.prodigal.faa", emit: predicted_proteins_faa
    path "{inputs.input_fasta.nameroot}.prodigal.ffn", emit: predicted_proteins_ffn
    path "{inputs.input_fasta.nameroot}.prodigal", emit: predicted_proteins_out

    script:
    def single_mode = params.prodigal_single_mode == false ? "" : "-p"
    """
    prodigal \
    -i ${input_fasta} \
    ${single_mode} \
    -a "${${input_fasta}.baseName}.prodigal.faa" \
    -d "${${input_fasta}.baseName}.prodigal.ffn" \
    -o "${${input_fasta}.baseName}.prodigal" \
    """

}


process PRODIGAL {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/prodigal"

    input:
    path input, stageAs: 'input'

    output:
    path "{inputs.identifier}.prodigal.ttl", emit: output

    script:
    """
    java -Xmx5g -jar /SAPP-2.0.jar -prodigal \
    -input ${input} \
    -output "${params.identifier}.prodigal.ttl" \
    """

}
