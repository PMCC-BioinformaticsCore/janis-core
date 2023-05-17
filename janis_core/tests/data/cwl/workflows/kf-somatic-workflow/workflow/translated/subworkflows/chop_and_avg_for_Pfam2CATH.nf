nextflow.enable.dsl=2

include { CHOP_AND_AVG_FROM_LIST } from './chop_and_avg_from_list'

workflow CHOP_AND_AVG_FOR_PFAM2_CATH {

    take:
    ch_crossmap_file
    ch_pdb_dir

    main:
    CHOP_AND_AVG_FROM_LIST(
        ch_crossmap_file.flatten().first(),
        ch_pdb_dir
    )

    emit:
    averaged_structs = CHOP_AND_AVG_FROM_LIST.out.avg_struct_per_fam

}
