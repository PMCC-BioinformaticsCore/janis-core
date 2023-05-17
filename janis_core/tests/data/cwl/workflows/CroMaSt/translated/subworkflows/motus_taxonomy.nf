nextflow.enable.dsl=2

include { CLEAN_CLASSIFICATION } from '../modules/clean_classification'
include { MOTUS_CLASSIFICATION } from '../modules/motus_classification'

workflow MOTUS_TAXONOMY {

    take:
    ch_reads

    main:
    CLEAN_CLASSIFICATION(
        MOTUS_CLASSIFICATION.out.motu_taxonomy
    )

    MOTUS_CLASSIFICATION(
        ch_reads
    )

    emit:
    motus = CLEAN_CLASSIFICATION.out.clean_annotations

}
