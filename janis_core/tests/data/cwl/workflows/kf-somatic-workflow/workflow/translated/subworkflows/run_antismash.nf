nextflow.enable.dsl=2

include { FIX_EMBL_AND_GBK } from '../modules/fix_embl_and_gbk'
include { FIX_GENECLUSTERS_TXT } from '../modules/fix_geneclusters_txt'
include { RUN_ANTISMASH } from '../modules/run_antismash'

workflow RUN_ANTISMASH {

    take:
    ch_accession
    ch_fasta_file
    ch_input_names_table

    main:
    FIX_EMBL_AND_GBK(
        RUN_ANTISMASH.out.embl_file,
        ch_input_names_table,
        ch_fasta_file,
        ch_fasta_file
    )

    FIX_GENECLUSTERS_TXT(
        RUN_ANTISMASH.out.geneclusters_txt,
        ch_fasta_file
    )

    RUN_ANTISMASH(
        ch_fasta_file,
        ch_accession,
        ch_fasta_file
    )

    emit:
    antismash_embl = FIX_EMBL_AND_GBK.out.fixed_embl
    antismash_gbk = FIX_EMBL_AND_GBK.out.fixed_gbk
    antismash_txt = FIX_GENECLUSTERS_TXT.out.fixed_txt

}
