nextflow.enable.dsl=2

include { EGGNOG } from './eggnog'
include { RUN_IPS } from './run_IPS'
include { RUN_HMMER } from './run_hmmer'
include { SPLIT_SEQS } from '../modules/split_seqs'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_cgc_predicted_proteins
    ch_egg_nog_data_dir
    ch_egg_nog_db
    ch_egg_nog_diamond_db
    ch_hmm_database
    ch_hmm_gathering_bit_score
    ch_hmm_omit_alignment
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_chunk_size_ips
    ch_chunk_size_eggnog
    ch_chunk_size_hmm
    ch_name_hmmer
    ch_name_ips

    main:
    EGGNOG(
        params.after_qc.functional_annotation_and_post_processing.functional_annotation.eggnog_cpu,
        ch_egg_nog_data_dir,
        ch_egg_nog_db,
        ch_egg_nog_diamond_db,
        SPLIT_SEQS.out.chunks,
        ch_cgc_predicted_proteins
    )

    RUN_IPS(
        ch_cgc_predicted_proteins,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_chunk_size_ips,
        ch_name_ips
    )

    RUN_HMMER(
        ch_cgc_predicted_proteins,
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_chunk_size_hmm,
        ch_name_hmmer,
        RUN_IPS.out.ips_result
    )

    SPLIT_SEQS(
        ch_cgc_predicted_proteins,
        ch_chunk_size_eggnog
    )

    emit:
    eggnog_annotations = EGGNOG.out.annotations
    eggnog_orthologs = EGGNOG.out.orthologs
    hmm_result = RUN_HMMER.out.hmm_result
    ips_result = RUN_IPS.out.ips_result

}


workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_cgc_predicted_proteins
    ch_hmm_database
    ch_hmm_gathering_bit_score
    ch_hmm_omit_alignment
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_chunk_size_ips
    ch_chunk_size_hmm
    ch_name_hmmer
    ch_name_ips

    main:
    RUN_IPS(
        ch_cgc_predicted_proteins,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_chunk_size_ips,
        ch_name_ips
    )

    RUN_IPS(
        ch_cgc_predicted_proteins,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases,
        ch_inter_pro_scan_output_format,
        ch_chunk_size_ips,
        ch_name_ips
    )

    RUN_HMMER(
        ch_cgc_predicted_proteins,
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_chunk_size_hmm,
        ch_name_hmmer,
        RUN_IPS.out.ips_result
    )

    RUN_HMMER(
        ch_cgc_predicted_proteins,
        ch_hmm_database,
        ch_hmm_gathering_bit_score,
        ch_hmm_omit_alignment,
        ch_chunk_size_hmm,
        ch_name_hmmer,
        RUN_IPS.out.ips_result
    )

    emit:
    hmm_result = RUN_HMMER.out.hmm_result
    ips_result = RUN_IPS.out.ips_result

}
