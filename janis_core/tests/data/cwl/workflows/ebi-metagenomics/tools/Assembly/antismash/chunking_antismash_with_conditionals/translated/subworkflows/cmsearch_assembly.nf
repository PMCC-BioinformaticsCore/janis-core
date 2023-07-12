nextflow.enable.dsl=2

include { CMSEARCH } from '../modules/cmsearch'
include { REMOVE_OVERLAPS } from '../modules/remove_overlaps'
include { RUN_CONCATENATE_MATCHES } from '../modules/run_concatenate_matches'

workflow CMSEARCH_ASSEMBLY {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences

    main:
    CMSEARCH(
        ch_query_sequences,
        ch_covariance_models.flatten().first()
    )

    REMOVE_OVERLAPS(
        RUN_CONCATENATE_MATCHES.out.result,
        ch_clan_info
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    emit:
    concatenate_matches = RUN_CONCATENATE_MATCHES.out.result
    deoverlapped_matches = REMOVE_OVERLAPS.out.deoverlapped_matches
    matches = CMSEARCH.out.matches

}


workflow CMSEARCH_ASSEMBLY {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences

    main:
    CMSEARCH(
        ch_query_sequences,
        ch_covariance_models.flatten().first()
    )

    CMSEARCH(
        ch_query_sequences,
        ch_covariance_models.flatten().first()
    )

    REMOVE_OVERLAPS(
        RUN_CONCATENATE_MATCHES.out.result,
        ch_clan_info
    )

    REMOVE_OVERLAPS(
        RUN_CONCATENATE_MATCHES.out.result,
        ch_clan_info
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    RUN_CONCATENATE_MATCHES(
        CMSEARCH.out.matches.toList(),
        ch_query_sequences
    )

    emit:
    concatenate_matches = RUN_CONCATENATE_MATCHES.out.result
    deoverlapped_matches = REMOVE_OVERLAPS.out.deoverlapped_matches
    matches = CMSEARCH.out.matches

}
