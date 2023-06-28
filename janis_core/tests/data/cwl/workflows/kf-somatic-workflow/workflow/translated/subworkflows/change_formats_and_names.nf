nextflow.enable.dsl=2

include { CREATE_CSV_GP } from '../modules/create_csv_gp'
include { CREATE_CSV_KP } from '../modules/create_csv_kp'

workflow CHANGE_FORMATS_AND_NAMES {

    take:
    ch_fasta
    ch_genome_properties_summary
    ch_kegg_summary

    main:
    CREATE_CSV_GP(
        ch_genome_properties_summary,
        ch_genome_properties_summary
    )

    CREATE_CSV_KP(
        ch_kegg_summary,
        ch_kegg_summary
    )

    emit:
    gp_summary_csv = CREATE_CSV_GP.out.csv_result
    kegg_summary_csv = CREATE_CSV_KP.out.csv_result

}
