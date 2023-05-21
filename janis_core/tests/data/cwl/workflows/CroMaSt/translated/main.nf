nextflow.enable.dsl=2

include { AFTER_QC } from './subworkflows/after_qc'
include { INDEX_STRELKA_BED } from './modules/index_strelka_bed'
include { PREPARE_REFERENCE } from './subworkflows/prepare_reference'
include { RUN_MANTA } from './subworkflows/run_manta'
include { NGTAX } from './modules/ngtax'
include { NGTAX_FILES_TO_FOLDER } from './modules/ngtax_files_to_folder'
include { NGTAX_TO_TSV_FASTA } from './modules/ngtax_to_tsv_fasta'
include { PHYLOSEQ_FILES_TO_FOLDER } from './modules/phyloseq_files_to_folder'
include { BAM_INDEX } from './modules/bam_index'
include { EXPRESSIONTOOL_BAM_INDEX } from './modules/expressiontool_bam_index'
include { PILON } from './modules/pilon'
include { READMAPPING_PILON } from './modules/readmapping_pilon'
include { CONVERSION } from './modules/conversion'
include { GZIP } from './modules/gzip'
include { INTERPROSCAN } from './modules/interproscan'
include { KOFAMSCAN } from './modules/kofamscan'
include { PRODIGAL } from './modules/prodigal'
include { COMPRESS } from './modules/compress'
include { TOHDT } from './modules/tohdt'
include { ADD_ISLAND_NAMES } from './modules/add_island_names'
include { ASSIGN_GENES } from './modules/assign_genes'
include { BED_TO_BIGBED } from './modules/bed_to_bigbed'
include { BED_TO_MACS } from './modules/bed_to_macs'
include { MAKE_GFF } from './modules/make_gff'
include { REDUCE_BED } from './modules/reduce_bed'
include { RENAME_PNG } from './modules/rename_png'
include { RUN_ROSE } from './modules/run_rose'
include { SORT_BED } from './modules/sort_bed'

ch_hg38_strelka_tbi     = Channel.fromPath( params.hg38_strelka_tbi ).ifEmpty( null )
ch_reference_dict       = Channel.fromPath( params.reference_dict ).ifEmpty( null )
ch_reference_fai        = Channel.fromPath( params.reference_fai ).ifEmpty( null )
ch_manta_cores          = Channel.of( params.manta_cores ).ifEmpty( null )
ch_manta_memory         = Channel.of( params.manta_memory ).ifEmpty( null )
ch_hg38_strelka_bed     = Channel.fromPath( params.hg38_strelka_bed )
ch_input_normal_aligned = Channel.fromPath( params.input_normal_aligned ).toList()
ch_input_tumor_aligned  = Channel.fromPath( params.input_tumor_aligned ).toList()
ch_reference_fasta      = Channel.fromPath( params.reference_fasta )
ch_vep_cache            = Channel.fromPath( params.vep_cache )


workflow  {

    INDEX_STRELKA_BED(
        ch_hg38_strelka_bed,
        ch_hg38_strelka_tbi
    )

    PREPARE_REFERENCE(
        ch_reference_fasta,
        ch_reference_dict,
        ch_reference_fai
    )

    RUN_MANTA(
        INDEX_STRELKA_BED.out.output,
        PREPARE_REFERENCE.out.indexed_fasta.map{ tuple -> tuple[0] },
        ch_input_normal_aligned,
        params.input_normal_name,
        ch_input_tumor_aligned,
        params.input_tumor_name,
        params.output_basename,
        PREPARE_REFERENCE.out.reference_dict,
        ch_vep_cache,
        ch_manta_cores,
        ch_manta_memory,
        params.select_vars_mode
    )


}


ch_forward_reads = Channel.fromPath( params.forward_reads )
ch_mapping_file  = Channel.fromPath( params.mapping_file )
ch_reverse_reads = Channel.fromPath( params.reverse_reads )


workflow  {

    NGTAX(
        ch_forward_reads,
        ch_mapping_file,
        ch_reverse_reads
    )


}


ch_metadata = Channel.fromPath( params.metadata ).ifEmpty( null )


workflow  {

    NGTAX()

    NGTAX()

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )


}


ch_embl = Channel.fromPath( params.embl )


workflow  {

    CONVERSION(
        ch_embl
    )

    GZIP(
        INTERPROSCAN.out.output
    )

    INTERPROSCAN(
        KOFAMSCAN.out.output
    )

    KOFAMSCAN(
        PRODIGAL.out.output
    )

    PRODIGAL(
        CONVERSION.out.output
    )


}


ch_input = Channel.fromPath( params.input )


workflow  {

    COMPRESS(
        TOHDT.out.output
    )

    TOHDT(
        ch_input
    )


}


ch_islands_control_file = Channel.fromPath( params.islands_control_file ).ifEmpty( null )
ch_annotation_file      = Channel.fromPath( params.annotation_file )
ch_bambai_pair          = Channel.fromPath( params.bambai_pair ).toList()
ch_chrom_length_file    = Channel.fromPath( params.chrom_length_file )
ch_islands_file         = Channel.fromPath( params.islands_file )


workflow  {

    ADD_ISLAND_NAMES(
        ASSIGN_GENES.out.result_file.toList(),
        ch_bambai_pair.map{ tuple -> tuple[0] }
    )

    ASSIGN_GENES(
        ch_annotation_file,
        BED_TO_MACS.out.output_file
    )

    BED_TO_BIGBED(
        ch_chrom_length_file,
        REDUCE_BED.out.output_file,
        ch_bambai_pair.map{ tuple -> tuple[0] }
    )

    BED_TO_MACS(
        SORT_BED.out.sorted_file
    )

    MAKE_GFF(
        ch_islands_file,
        ch_islands_control_file
    )

    REDUCE_BED(
        SORT_BED.out.sorted_file
    )

    RENAME_PNG(
        RUN_ROSE.out.plot_points_pic,
        ch_bambai_pair.map{ tuple -> tuple[0] }
    )

    RUN_ROSE(
        ch_bambai_pair,
        ch_annotation_file,
        MAKE_GFF.out.gff_file
    )

    SORT_BED(
        RUN_ROSE.out.gateway_super_enhancers_bed
    )


}
