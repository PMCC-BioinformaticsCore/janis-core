nextflow.enable.dsl=2

include { AVG_AVERAGED_STRUCTURES } from '../modules/avg_averaged_structures'
include { AVG_UNP_DOMAINS } from './avg_unp_domains'
include { COPY_AVG_DOM } from '../modules/copy_avg_dom'
include { PER_UNP_DOM_INSTANCE } from '../modules/per_unp_dom_instance'

workflow CHOP_AND_AVG_FOR_CORE {

    take:
    ch_core_list
    ch_pdb_dir

    main:
    AVG_AVERAGED_STRUCTURES(
        PER_UNP_DOM_INSTANCE.out.family_name,
        COPY_AVG_DOM.out.dir_unp_dom
    )

    AVG_UNP_DOMAINS(
        PER_UNP_DOM_INSTANCE.out.dom_per_fam.flatten().first(),
        ch_pdb_dir
    )

    COPY_AVG_DOM(
        AVG_UNP_DOMAINS.out.avg_unp_dom_structs.toList()
    )

    PER_UNP_DOM_INSTANCE(
        ch_core_list
    )

    emit:
    averaged_structs = AVG_AVERAGED_STRUCTURES.out.avg_structs

}
