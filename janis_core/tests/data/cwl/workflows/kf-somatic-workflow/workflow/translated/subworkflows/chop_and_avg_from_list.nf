nextflow.enable.dsl=2

include { AVG_AVERAGED_UNP_DOMAINS } from '../modules/avg_averaged_unp_domains'
include { AVG_UNP_DOMAINS } from './avg_unp_domains'
include { COPY_AVG_DOM } from '../modules/copy_avg_dom'
include { PER_DOM_INSTANCE } from '../modules/per_dom_instance'

workflow CHOP_AND_AVG_FROM_LIST {

    take:
    ch_in_file
    ch_pdb_storage

    main:
    AVG_AVERAGED_UNP_DOMAINS(
        PER_DOM_INSTANCE.out.family_name,
        COPY_AVG_DOM.out.dir_unp_dom
    )

    AVG_UNP_DOMAINS(
        PER_DOM_INSTANCE.out.dom_per_fam.flatten().first(),
        ch_pdb_storage
    )

    COPY_AVG_DOM(
        AVG_UNP_DOMAINS.out.avg_structures.toList()
    )

    PER_DOM_INSTANCE(
        ch_in_file
    )

    emit:
    avg_struct_per_fam = AVG_AVERAGED_UNP_DOMAINS.out.avg_structs

}
