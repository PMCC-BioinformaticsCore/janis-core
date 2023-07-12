nextflow.enable.dsl=2

include { ADD_DOMAIN_POSITIONS } from '../modules/add_domain_positions'
include { COLLECT_LOST_INSTANCES } from '../modules/collect_lost_instances'
include { FILTER_PFAM_STRUCTURES } from '../modules/filter_pfam_structures'
include { RESMAPPING_PFAM_STRUCTS } from './resmapping_pfam_structs'

workflow PFAM_DOMAIN_INSTANCES {

    take:
    ch_family_idsfile
    ch_lost_merged
    ch_min_dom_size
    ch_resmapped_file
    ch_siftsdir

    main:
    ADD_DOMAIN_POSITIONS(
        RESMAPPING_PFAM_STRUCTS.out.resmapped_pfam.toList(),
        ch_resmapped_file
    )

    COLLECT_LOST_INSTANCES(
        RESMAPPING_PFAM_STRUCTS.out.lost_insta_pfam.toList(),
        FILTER_PFAM_STRUCTURES.out.pfam_obs,
        ch_lost_merged
    )

    FILTER_PFAM_STRUCTURES(
        ch_family_idsfile,
        ch_min_dom_size
    )

    RESMAPPING_PFAM_STRUCTS(
        FILTER_PFAM_STRUCTURES.out.splitted_pfam_sep.flatten().first(),
        ch_siftsdir
    )

    emit:
    pfam_domain_posi_file = ADD_DOMAIN_POSITIONS.out.resmapped_domains
    pfam_total_lost_structures = COLLECT_LOST_INSTANCES.out.lost_domain_list

}
