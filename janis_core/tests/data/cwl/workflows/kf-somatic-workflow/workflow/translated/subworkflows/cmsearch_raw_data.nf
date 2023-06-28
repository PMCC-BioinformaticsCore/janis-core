nextflow.enable.dsl=2

include { CAT_MODELS } from '../modules/cat_models'
include { CMSEARCH } from '../modules/cmsearch'
include { REMOVE_OVERLAPS } from '../modules/remove_overlaps'
include { RUN_CONCATENATE_DEOVERLAPPED_MATCHES } from '../modules/run_concatenate_deoverlapped_matches'
include { RUN_CONCATENATE_MATCHES } from '../modules/run_concatenate_matches'
include { SPLIT_FASTA } from '../modules/split_fasta'

workflow CMSEARCH_RAW_DATA {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences

    main:
    CAT_MODELS(
        ch_covariance_models
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    emit:
    concatenate_matches = RUN_CONCATENATE_MATCHES.out.result
    deoverlapped_matches = RUN_CONCATENATE_DEOVERLAPPED_MATCHES.out.result

}


workflow CMSEARCH_RAW_DATA {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences

    main:
    CAT_MODELS(
        ch_covariance_models
    )

    CAT_MODELS(
        ch_covariance_models
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    emit:
    concatenate_matches = RUN_CONCATENATE_MATCHES.out.result
    deoverlapped_matches = RUN_CONCATENATE_DEOVERLAPPED_MATCHES.out.result

}


workflow CMSEARCH_RAW_DATA {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences

    main:
    CAT_MODELS(
        ch_covariance_models
    )

    CAT_MODELS(
        ch_covariance_models
    )

    CAT_MODELS(
        ch_covariance_models
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    CMSEARCH(
        SPLIT_FASTA.out.chunks.flatten().first(),
        CAT_MODELS.out.result
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    REMOVE_OVERLAPS(
        CMSEARCH.out.matches,
        ch_clan_info
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_DEOVERLAPPED_MATCHES(
        REMOVE_OVERLAPS.out.deoverlapped_matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    SPLIT_FASTA(
        ch_query_sequences
    )

    emit:
    concatenate_matches = RUN_CONCATENATE_MATCHES.out.result
    deoverlapped_matches = RUN_CONCATENATE_DEOVERLAPPED_MATCHES.out.result

}
