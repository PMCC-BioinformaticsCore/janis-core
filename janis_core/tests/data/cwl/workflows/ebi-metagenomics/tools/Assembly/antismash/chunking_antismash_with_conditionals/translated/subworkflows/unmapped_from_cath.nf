nextflow.enable.dsl=2

include { AVG_UNP_DOMAINS } from './avg_unp_domains'
include { CHECK_THRESHOLD_STEP } from '../modules/check_threshold_step'
include { COPY_AVG_DOM } from '../modules/copy_avg_dom'
include { PAIRWISE_ALIGN_AVG_STRUCTS } from '../modules/pairwise_align_avg_structs'
include { PER_UNP_DOM_INSTANCE } from '../modules/per_unp_dom_instance'

workflow UNMAPPED_FROM_CATH {

    take:
    ch_alignment_score
    ch_core_struct
    ch_failed_name
    ch_iteration
    ch_passed_name
    ch_pdb_dir
    ch_score_threshold
    ch_unmapped_list

    main:
    AVG_UNP_DOMAINS(
        PER_UNP_DOM_INSTANCE.out.dom_per_fam.flatten().first(),
        ch_pdb_dir
    )

    CHECK_THRESHOLD_STEP(
        ch_unmapped_list,
        PAIRWISE_ALIGN_AVG_STRUCTS.out.alignment_out,
        ch_failed_name,
        ch_passed_name,
        ch_alignment_score,
        ch_score_threshold
    )

    COPY_AVG_DOM(
        AVG_UNP_DOMAINS.out.avg_unp_dom_structs.toList()
    )

    PAIRWISE_ALIGN_AVG_STRUCTS(
        ch_core_struct,
        COPY_AVG_DOM.out.dir_unp_dom
    )

    PER_UNP_DOM_INSTANCE(
        ch_unmapped_list
    )

    emit:
    domain_like_list = CHECK_THRESHOLD_STEP.out.passed_structs_list
    failed_domains_list = CHECK_THRESHOLD_STEP.out.failed_structs_list
    unmapped_aligned_results = PAIRWISE_ALIGN_AVG_STRUCTS.out.alignment_out

}
