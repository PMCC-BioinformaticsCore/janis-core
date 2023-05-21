nextflow.enable.dsl=2

include { RESMAPPING_FOR_PFAM_UP2_PDB } from '../modules/resmapping_for_Pfam_UP2PDB'

workflow RESMAPPING_PFAM_STRUCTS {

    take:
    ch_flt_files
    ch_sifts

    main:
    RESMAPPING_FOR_PFAM_UP2_PDB(
        ch_flt_files,
        ch_sifts
    )

    emit:
    lost_insta_pfam = RESMAPPING_FOR_PFAM_UP2_PDB.out.pfam_lost
    resmapped_pfam = RESMAPPING_FOR_PFAM_UP2_PDB.out.pfam_resmapped

}
