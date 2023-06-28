nextflow.enable.dsl=2

include { RESMAPPING_FOR_CATH_PDB2_UP } from '../modules/resmapping_for_CATH_PDB2UP'

workflow RESMAPPING_CATH_STRUCTS {

    take:
    ch_flt_files
    ch_sifts

    main:
    RESMAPPING_FOR_CATH_PDB2_UP(
        ch_flt_files,
        ch_sifts
    )

    emit:
    lost_insta_cath = RESMAPPING_FOR_CATH_PDB2_UP.out.cath_lost
    resmapped_cath = RESMAPPING_FOR_CATH_PDB2_UP.out.cath_resmapped

}
