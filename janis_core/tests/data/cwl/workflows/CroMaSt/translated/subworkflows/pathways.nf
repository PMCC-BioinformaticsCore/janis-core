nextflow.enable.dsl=2

include { KEGG_PATHWAYS } from '../modules/kegg_pathways'
include { PARSING_HMMSCAN } from '../modules/parsing_hmmscan'
include { TAB_MODIFICATION } from '../modules/tab_modification'

workflow PATHWAYS {

    take:
    ch_filtered_fasta
    ch_graphs
    ch_input_table_hmmscan
    ch_outputname
    ch_pathways_classes
    ch_pathways_names

    main:
    KEGG_PATHWAYS(
        PARSING_HMMSCAN.out.output_table,
        ch_graphs,
        ch_outputname,
        ch_pathways_classes,
        ch_pathways_names
    )

    PARSING_HMMSCAN(
        ch_filtered_fasta,
        TAB_MODIFICATION.out.output_with_tabs
    )

    TAB_MODIFICATION(
        ch_input_table_hmmscan
    )

    emit:
    kegg_contigs_summary = KEGG_PATHWAYS.out.summary_contigs
    kegg_pathways_summary = KEGG_PATHWAYS.out.summary_pathways
    modification_out = TAB_MODIFICATION.out.output_with_tabs
    parsing_hmmscan_out = PARSING_HMMSCAN.out.output_table

}
