nextflow.enable.dsl=2

include { ADD_ISLAND_NAMES_TOOL as ADD_ISLAND_NAMES } from './modules/add_island_names_tool'
include { ASSIGN_GENES_TOOL as ASSIGN_GENES } from './modules/assign_genes_tool'
include { BED_TO_BIGBED_TOOL as BED_TO_BIGBED } from './modules/bed_to_bigbed_tool'
include { BED_TO_MACS_TOOL as BED_TO_MACS } from './modules/bed_to_macs_tool'
include { MAKE_GFF_TOOL as MAKE_GFF } from './modules/make_gff_tool'
include { REDUCE_BED_TOOL as REDUCE_BED } from './modules/reduce_bed_tool'
include { RENAME_PNG_TOOL as RENAME_PNG } from './modules/rename_png_tool'
include { RUN_ROSE_TOOL as RUN_ROSE } from './modules/run_rose_tool'
include { SORT_BED_TOOL as SORT_BED } from './modules/sort_bed_tool'


// data which will be passed as variables
annotation_file       = file( params.annotation_file )
bambai_pair           = params.bambai_pair.collect{ file(it) }
chrom_length_file     = file( params.chrom_length_file )
islands_control_file  = file( params.islands_control_file )
islands_file          = file( params.islands_file )


workflow {

    ADD_ISLAND_NAMES(
        ASSIGN_GENES.out.result_file.toList(),  // input_file
        bambai_pair.map{ tuple -> tuple[0] }    // param
    )

    ASSIGN_GENES(
        annotation_file,              // annotation_filename
        BED_TO_MACS.out.output_file,  // input_filename
        params.promoter_bp            // promoter_bp
    )

    BED_TO_BIGBED(
        chrom_length_file,                    // chrom_length_file
        REDUCE_BED.out.output_file,           // input_bed
        "bed4",                               // bed_type
        bambai_pair.map{ tuple -> tuple[0] }  // output_filename
    )

    BED_TO_MACS(
        SORT_BED.out.sorted_file  // input_file
    )

    MAKE_GFF(
        islands_file,         // islands_file
        islands_control_file  // islands_control_file
    )

    REDUCE_BED(
        SORT_BED.out.sorted_file  // input_file
    )

    RENAME_PNG(
        RUN_ROSE.out.plot_points_pic,         // source_file
        bambai_pair.map{ tuple -> tuple[0] }  // target_filename
    )

    RUN_ROSE(
        bambai_pair,             // bam_file
        annotation_file,         // annotation_file
        MAKE_GFF.out.gff_file,   // binding_sites_file
        params.stitch_distance,  // stitch_distance
        params.tss_distance      // tss_distance
    )

    SORT_BED(
        RUN_ROSE.out.gateway_super_enhancers_bed,  // unsorted_file
        ['1,1', '2,2n', '3,3n']                    // key
    )


}
