nextflow.enable.dsl=2

include { ADD_DOMAIN_POSITIONS } from '../modules/add_domain_positions'
include { COLLECT_LOST_INSTANCES } from '../modules/collect_lost_instances'
include { FILTER_CATH_STRUCTURES } from '../modules/filter_cath_structures'
include { RESMAPPING_CATH_STRUCTS } from './resmapping_cath_structs'

workflow CATH_DOMAIN_INSTANCES {

    take:
    ch_family_idsfile
    ch_lost_merged
    ch_min_dom_size
    ch_resmapped_file
    ch_siftsdir

    main:
    ADD_DOMAIN_POSITIONS(
        RESMAPPING_CATH_STRUCTS.out.resmapped_cath.toList(),
        ch_resmapped_file
    )

    COLLECT_LOST_INSTANCES(
        RESMAPPING_CATH_STRUCTS.out.lost_insta_cath.toList(),
        FILTER_CATH_STRUCTURES.out.cath_obs,
        ch_lost_merged
    )

    FILTER_CATH_STRUCTURES(
        ch_family_idsfile,
        ch_min_dom_size
    )

    RESMAPPING_CATH_STRUCTS(
        FILTER_CATH_STRUCTURES.out.splitted_cath_sep.flatten().first(),
        ch_siftsdir
    )

    emit:
    cath_domain_posi_file = ADD_DOMAIN_POSITIONS.out.resmapped_domains
    cath_total_lost_structures = COLLECT_LOST_INSTANCES.out.lost_domain_list

}
