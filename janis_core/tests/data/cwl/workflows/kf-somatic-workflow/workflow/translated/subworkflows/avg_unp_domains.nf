nextflow.enable.dsl=2

include { AVG_CHOPPED_STRUCTS_UNP_DOMAINS } from '../modules/avg_chopped_structs_unp_domains'
include { CHOP_STRUCTS } from '../modules/chop_structs'

workflow AVG_UNP_DOMAINS {

    take:
    ch_domfiles
    ch_pdb_d

    main:
    AVG_CHOPPED_STRUCTS_UNP_DOMAINS(
        CHOP_STRUCTS.out.family_name,
        CHOP_STRUCTS.out.split_structs_dir
    )

    CHOP_STRUCTS(
        ch_domfiles,
        ch_pdb_d
    )

    emit:
    avg_unp_dom_structs = AVG_CHOPPED_STRUCTS_UNP_DOMAINS.out.avg_structs

}
