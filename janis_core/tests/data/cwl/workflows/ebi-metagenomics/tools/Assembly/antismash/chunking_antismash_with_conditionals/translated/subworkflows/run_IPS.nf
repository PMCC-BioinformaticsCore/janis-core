nextflow.enable.dsl=2

include { COMBINE_IPS } from '../modules/combine_ips'
include { INTERPROSCAN } from '../modules/interproscan'
include { SPLIT_SEQS } from '../modules/split_seqs'

workflow RUN_IPS {

    take:
    ch_cgc_predicted_proteins
    ch_inter_pro_scan_applications
    ch_inter_pro_scan_databases
    ch_inter_pro_scan_output_format
    ch_chunk_size
    ch_name_ips

    main:
    COMBINE_IPS(
        INTERPROSCAN.out.i5Annotations.toList(),
        ch_cgc_predicted_proteins,
        ch_name_ips
    )

    INTERPROSCAN(
        SPLIT_SEQS.out.chunks.flatten().first(),
        ch_inter_pro_scan_output_format,
        ch_inter_pro_scan_applications,
        ch_inter_pro_scan_databases
    )

    SPLIT_SEQS(
        ch_cgc_predicted_proteins,
        ch_chunk_size
    )

    emit:
    ips_result = COMBINE_IPS.out.result

}
