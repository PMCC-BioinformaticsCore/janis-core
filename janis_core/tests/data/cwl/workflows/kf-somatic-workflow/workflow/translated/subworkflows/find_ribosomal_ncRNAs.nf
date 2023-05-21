nextflow.enable.dsl=2

include { CMSEARCH_ASSEMBLY } from './cmsearch_assembly'
include { CMSEARCH_RAW_DATA } from './cmsearch_raw_data'

workflow FIND_RIBOSOMAL_NC_RNAS {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences
    ch_type

    main:
    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    emit:
    concatenate_matches = [CMSEARCH_ASSEMBLY.out.concatenate_matches, CMSEARCH_RAW_DATA.out.concatenate_matches]
    deoverlapped_matches = [CMSEARCH_ASSEMBLY.out.deoverlapped_matches, CMSEARCH_RAW_DATA.out.deoverlapped_matches]

}


workflow FIND_RIBOSOMAL_NC_RNAS {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences
    ch_type

    main:
    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    emit:
    concatenate_matches = [CMSEARCH_ASSEMBLY.out.concatenate_matches, CMSEARCH_RAW_DATA.out.concatenate_matches]
    deoverlapped_matches = [CMSEARCH_ASSEMBLY.out.deoverlapped_matches, CMSEARCH_RAW_DATA.out.deoverlapped_matches]

}


workflow FIND_RIBOSOMAL_NC_RNAS {

    take:
    ch_clan_info
    ch_covariance_models
    ch_query_sequences
    ch_type

    main:
    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_ASSEMBLY(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    CMSEARCH_RAW_DATA(
        ch_clan_info,
        ch_covariance_models,
        ch_query_sequences
    )

    emit:
    concatenate_matches = [CMSEARCH_ASSEMBLY.out.concatenate_matches, CMSEARCH_RAW_DATA.out.concatenate_matches]
    deoverlapped_matches = [CMSEARCH_ASSEMBLY.out.deoverlapped_matches, CMSEARCH_RAW_DATA.out.deoverlapped_matches]

}
