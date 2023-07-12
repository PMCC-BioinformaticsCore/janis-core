nextflow.enable.dsl=2

process PICARD_CREATE_SEQUENCE_DICTIONARY {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.7.0R"
    publishDir "${params.outdir}/prepare_reference/picard_create_sequence_dictionary"

    input:
    path input_fasta, stageAs: 'input_fasta'
    path input_dict, stageAs: 'input_dict'

    output:
    path "*.dict", emit: dict

    script:
    """
    \
    <js>inputs.input_dict ? 'echo java -jar /gatk-package-4.1.7.0-local.jar' : 'java -jar /gatk-package-4.1.7.0-local.jar' </js> \
    CreateSequenceDictionary \
    -R ${input_fasta} \
    """

}
