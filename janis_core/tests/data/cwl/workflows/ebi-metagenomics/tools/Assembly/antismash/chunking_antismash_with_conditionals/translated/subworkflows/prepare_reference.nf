nextflow.enable.dsl=2

include { BUNDLE_SECONDARIES } from '../modules/bundle_secondaries'
include { BWA_INDEX } from '../modules/bwa_index'
include { PICARD_CREATE_SEQUENCE_DICTIONARY } from '../modules/picard_create_sequence_dictionary'
include { SAMTOOLS_FAIDX } from '../modules/samtools_faidx'

workflow PREPARE_REFERENCE {

    take:
    ch_input_fasta
    ch_input_dict
    ch_input_fai

    main:
    BUNDLE_SECONDARIES(
        ch_input_fasta,
        SAMTOOLS_FAIDX.out.fai.toList()
    )

    BWA_INDEX(
        ch_input_fasta
    )

    PICARD_CREATE_SEQUENCE_DICTIONARY(
        ch_input_fasta,
        ch_input_dict
    )

    SAMTOOLS_FAIDX(
        ch_input_fasta,
        ch_input_fai
    )

    emit:
    indexed_fasta = BUNDLE_SECONDARIES.out.output
    reference_dict = PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict

}
