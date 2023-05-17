nextflow.enable.dsl=2

include { TP_CUT_TOOL as TP_CUT_TOOL1 } from './modules/tp_cut_tool'
include { MERGE_COLS1 } from './modules/merge_cols1'
include { TP_REPLACE_IN_LINE } from './modules/tp_replace_in_line'
include { TP_CUT_TOOL as TP_CUT_TOOL2 } from './modules/tp_cut_tool'
include { ANNOTATEMYIDS } from './modules/annotatemyids'
include { LIMMA_VOOM } from './modules/limma_voom'


// data which will be passed as channels
ch_in_sampleinfo                 = Channel.fromPath( params.in_sampleinfo )
ch_in_seqdata                    = Channel.fromPath( params.in_seqdata )
ch_limma_voom_limma_voom_script  = Channel.fromPath( params.limma_voom_limma_voom_script )


workflow {

    TP_CUT_TOOL1(
        ch_in_seqdata  // inputFile
    )

    MERGE_COLS1(
        ch_in_sampleinfo  // input1
    )

    TP_REPLACE_IN_LINE(
        TP_CUT_TOOL1.out.outputFile  // infile
    )

    TP_CUT_TOOL2(
        MERGE_COLS1.out.out_file12  // inputFile
    )

    ANNOTATEMYIDS(
        TP_REPLACE_IN_LINE.out.outfile  // unknown1
    )

    LIMMA_VOOM(
        ch_limma_voom_limma_voom_script,  // limma_voom_script
        ANNOTATEMYIDS.out.out_tab,        // option_a
        TP_CUT_TOOL2.out.outputFile,      // option_f
        TP_REPLACE_IN_LINE.out.outfile    // option_m
    )


}
