nextflow.enable.dsl=2

include { ADD_CROSSMAPPED_TO_RESMAPPED } from './modules/add_crossmapped_to_resmapped'
include { ALIGN_AVG_STRUCTS_PAIRWISE } from './modules/align_avg_structs_pairwise'
include { CATH_DOMAIN_INSTANCES } from './subworkflows/cath_domain_instances'
include { CHECK_ALIGNMENT_SCORES } from './modules/check_alignment_scores'
include { CHOP_AND_AVG_FOR_CATH2PFAM } from './subworkflows/chop_and_avg_for_CATH2Pfam'
include { CHOP_AND_AVG_FOR_PFAM2_CATH } from './subworkflows/chop_and_avg_for_Pfam2CATH'
include { CHOP_AND_AVG_FOR_CORE } from './subworkflows/chop_and_avg_for_core'
include { COMPARE_INSTANCES_CATH_PFAM } from './modules/compare_instances_CATH_Pfam'
include { CREATE_NEW_PARAMETERS } from './modules/create_new_parameters'
include { CROSSMAPPING_CATH2PFAM } from './modules/crossmapping_CATH2Pfam'
include { CROSSMAPPING_PFAM2_CATH } from './modules/crossmapping_Pfam2CATH'
include { FORMAT_CORE_LIST } from './modules/format_core_list'
include { GATHER_DOMAIN_LIKE } from './modules/gather_domain_like'
include { GATHER_FAILED_DOMAINS } from './modules/gather_failed_domains'
include { GET_FAMILY_IDS } from './modules/get_family_ids'
include { PFAM_DOMAIN_INSTANCES } from './subworkflows/pfam_domain_instances'
include { UNMAPPED_FROM_CATH } from './subworkflows/unmapped_from_cath'
include { UNMAPPED_FROM_PFAM } from './subworkflows/unmapped_from_pfam'
include { AFTER_QC } from './subworkflows/after_qc'
include { BEFORE_QC } from './subworkflows/before_qc'
include { TOUCH_FILE_FLAG } from './modules/touch_file_flag'
include { TOUCH_NO_CDS_FLAG } from './modules/touch_no_cds_flag'
include { INDEX_STRELKA_BED } from './modules/index_strelka_bed'
include { PREPARE_REFERENCE } from './subworkflows/prepare_reference'
include { RUN_MANTA } from './subworkflows/run_manta'
include { NGTAX } from './modules/ngtax'
include { CARVEME } from './modules/carveme'
include { CARVEME_FILES_TO_FOLDER } from './modules/carveme_files_to_folder'
include { COMPRESS_CARVEME } from './modules/compress_carveme'
include { COMPRESS_PRODIGAL } from './modules/compress_prodigal'
include { GEMSTATS } from './modules/gemstats'
include { MEMOTE_FILES_TO_FOLDER } from './modules/memote_files_to_folder'
include { MEMOTE_REPORT_SNAPSHOT } from './modules/memote_report_snapshot'
include { MEMOTE_RUN } from './modules/memote_run'
include { PRODIGAL } from './modules/prodigal'
include { PRODIGAL_FILES_TO_FOLDER } from './modules/prodigal_files_to_folder'
include { SMETANA } from './modules/smetana'
include { NGTAX_FILES_TO_FOLDER } from './modules/ngtax_files_to_folder'
include { NGTAX_TO_TSV_FASTA } from './modules/ngtax_to_tsv_fasta'
include { PHYLOSEQ_FILES_TO_FOLDER } from './modules/phyloseq_files_to_folder'
include { FASTQC } from './modules/fastqc'
include { FASTQC_FILES_TO_FOLDER } from './modules/fastqc_files_to_folder'
include { READS_TO_FOLDER } from './modules/reads_to_folder'
include { FOLDER_COMPRESSION } from './modules/folder_compression'
include { PICRUST2 } from './modules/picrust2'
include { PICRUST_FILES_TO_FOLDER } from './modules/picrust_files_to_folder'
include { BAM_INDEX } from './modules/bam_index'
include { EXPRESSIONTOOL_BAM_INDEX } from './modules/expressiontool_bam_index'
include { PILON } from './modules/pilon'
include { READMAPPING_PILON } from './modules/readmapping_pilon'
include { SAM_TO_SORTED_BAM } from './modules/sam_to_sorted_bam'
include { VCF_COMPRESS } from './modules/vcf_compress'
include { CONVERSION } from './modules/conversion'
include { GZIP } from './modules/gzip'
include { INTERPROSCAN } from './modules/interproscan'
include { KOFAMSCAN } from './modules/kofamscan'
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

ch_core_domain_struct     = Channel.fromPath( params.core_domain_struct )
ch_domain_like            = Channel.fromPath( params.domain_like )
ch_failed_domain          = Channel.fromPath( params.failed_domain )
ch_filename               = Channel.fromPath( params.filename )
ch_paramfile              = Channel.fromPath( params.paramfile )
ch_prev_cross_mapped_cath = Channel.fromPath( params.prev_cross_mapped_cath )
ch_prev_cross_mapped_pfam = Channel.fromPath( params.prev_cross_mapped_pfam )
ch_true_domain_file       = Channel.fromPath( params.true_domain_file )


workflow  {

    ADD_CROSSMAPPED_TO_RESMAPPED(
        ch_prev_cross_mapped_cath,
        CATH_DOMAIN_INSTANCES.out.cath_domain_posi_file,
        ch_prev_cross_mapped_pfam,
        PFAM_DOMAIN_INSTANCES.out.pfam_domain_posi_file
    )

    ALIGN_AVG_STRUCTS_PAIRWISE(
        CHOP_AND_AVG_FOR_CATH2PFAM.out.averaged_structs.flatten().first(),
        CHOP_AND_AVG_FOR_CORE.out.averaged_structs,
        CHOP_AND_AVG_FOR_PFAM2_CATH.out.averaged_structs.flatten().first()
    )

    CATH_DOMAIN_INSTANCES(
        GET_FAMILY_IDS.out.family_ids,
        params.cath_lost,
        params.min_domain_length,
        params.cath_resmap,
        params.sifts_dir
    )

    CHECK_ALIGNMENT_SCORES(
        ALIGN_AVG_STRUCTS_PAIRWISE.out.alignment_out,
        GET_FAMILY_IDS.out.family_ids,
        CROSSMAPPING_CATH2PFAM.out.allcrossmap_cath,
        CROSSMAPPING_PFAM2_CATH.out.allcrossmap_pfam
    )

    CHOP_AND_AVG_FOR_CATH2PFAM(
        CROSSMAPPING_CATH2PFAM.out.cath_crossmapped,
        params.pdb_dir
    )

    CHOP_AND_AVG_FOR_PFAM2_CATH(
        CROSSMAPPING_PFAM2_CATH.out.pfam_crossmapped,
        params.pdb_dir
    )

    CHOP_AND_AVG_FOR_CORE(
        FORMAT_CORE_LIST.out.coredomains_list,
        params.pdb_dir
    )

    COMPARE_INSTANCES_CATH_PFAM(
        CATH_DOMAIN_INSTANCES.out.cath_domain_posi_file,
        PFAM_DOMAIN_INSTANCES.out.pfam_domain_posi_file,
        ch_true_domain_file
    )

    CREATE_NEW_PARAMETERS(
        CHOP_AND_AVG_FOR_CORE.out.averaged_structs,
        CHECK_ALIGNMENT_SCORES.out.cath_crossmap_passed,
        CHECK_ALIGNMENT_SCORES.out.pfam_crossmap_passed,
        GATHER_DOMAIN_LIKE.out.unmapped_list,
        GATHER_FAILED_DOMAINS.out.unmapped_list,
        GET_FAMILY_IDS.out.family_ids,
        ch_paramfile,
        COMPARE_INSTANCES_CATH_PFAM.out.common_domains
    )

    CROSSMAPPING_CATH2PFAM(
        COMPARE_INSTANCES_CATH_PFAM.out.cath_unique
    )

    CROSSMAPPING_PFAM2_CATH(
        COMPARE_INSTANCES_CATH_PFAM.out.pfam_unique
    )

    FORMAT_CORE_LIST(
        COMPARE_INSTANCES_CATH_PFAM.out.common_domains
    )

    GATHER_DOMAIN_LIKE(
        ch_domain_like,
        UNMAPPED_FROM_CATH.out.domain_like_list,
        UNMAPPED_FROM_PFAM.out.domain_like_list
    )

    GATHER_FAILED_DOMAINS(
        ch_failed_domain,
        CHECK_ALIGNMENT_SCORES.out.cath_crossmap_failed,
        UNMAPPED_FROM_CATH.out.failed_domains_list,
        CHECK_ALIGNMENT_SCORES.out.pfam_crossmap_failed,
        UNMAPPED_FROM_PFAM.out.failed_domains_list
    )

    GET_FAMILY_IDS(
        ch_filename
    )

    PFAM_DOMAIN_INSTANCES(
        GET_FAMILY_IDS.out.family_ids,
        params.pfam_lost,
        params.min_domain_length,
        params.pfam_resmap,
        params.sifts_dir
    )

    UNMAPPED_FROM_CATH(
        params.alignment_score,
        CHOP_AND_AVG_FOR_CORE.out.averaged_structs,
        params.unmap_cath_fail,
        params.iteration,
        params.unmap_cath_pass,
        params.pdb_dir,
        params.score_threshold,
        CROSSMAPPING_CATH2PFAM.out.cath_unmapped
    )

    UNMAPPED_FROM_PFAM(
        params.alignment_score,
        CHOP_AND_AVG_FOR_CORE.out.averaged_structs,
        params.unmap_pfam_fail,
        params.iteration,
        params.unmap_pfam_pass,
        params.pdb_dir,
        params.score_threshold,
        CROSSMAPPING_PFAM2_CATH.out.pfam_unmapped
    )


}


ch_forward_reads = Channel.fromPath( params.forward_reads ).ifEmpty( null )
ch_reverse_reads = Channel.fromPath( params.reverse_reads ).ifEmpty( null )
ch_single_reads  = Channel.fromPath( params.single_reads ).ifEmpty( null )
ch_itsonedb      = Channel.fromPath( params.itsonedb ).toList()
ch_lsu_db        = Channel.fromPath( params.lsu_db ).toList()
ch_ssu_db        = Channel.fromPath( params.ssu_db ).toList()
ch_unite_db      = Channel.fromPath( params.unite_db ).toList()


workflow  {

    AFTER_QC(
        params.5.8s_pattern,
        params.5s_pattern,
        BEFORE_QC.out.filtered_fasta,
        ch_itsonedb,
        params.itsonedb_label,
        params.itsonedb_otu_file,
        params.itsonedb_tax,
        ch_lsu_db,
        params.lsu_label,
        params.lsu_otus,
        params.lsu_tax,
        params.rfam_model_clans,
        params.rfam_models,
        ch_ssu_db,
        params.ssu_label,
        params.ssu_otus,
        params.ssu_tax,
        ch_unite_db,
        params.unite_label,
        params.unite_otu_file,
        params.unite_tax
    )

    BEFORE_QC(
        params.qc_min_length,
        params.stats_file_name,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    TOUCH_FILE_FLAG()


}


ch_forward_reads      = Channel.fromPath( params.forward_reads ).ifEmpty( null )
ch_reverse_reads      = Channel.fromPath( params.reverse_reads ).ifEmpty( null )
ch_single_reads       = Channel.fromPath( params.single_reads ).ifEmpty( null )
ch_egg_nog_data_dir   = Channel.of( params.egg_nog_data_dir ).ifEmpty( null )
ch_egg_nog_db         = Channel.of( params.egg_nog_db ).ifEmpty( null )
ch_egg_nog_diamond_db = Channel.of( params.egg_nog_diamond_db ).ifEmpty( null )
ch_lsu_db             = Channel.fromPath( params.lsu_db ).toList()
ch_ssu_db             = Channel.fromPath( params.ssu_db ).toList()


workflow  {

    AFTER_QC(
        params.5.8s_pattern,
        params.5s_pattern,
        params.cgc_config,
        params.cgc_postfixes,
        params.hmm_gathering_bit_score,
        params.hmm_name_database,
        params.hmm_omit_alignment,
        params.inter_pro_scan_applications,
        params.inter_pro_scan_databases,
        params.inter_pro_scan_output_format,
        params.cgc_chunk_size,
        BEFORE_QC.out.filtered_fasta,
        params.func_ann_names_hmmer,
        params.func_ann_names_ips,
        params.go_config,
        params.hmmsearch_header,
        params.ips_header,
        params.ko_file,
        ch_lsu_db,
        params.lsu_label,
        params.lsu_otus,
        params.lsu_tax,
        BEFORE_QC.out.motus_input,
        params.other_ncrna_models,
        params.protein_chunk_size_ips,
        params.protein_chunk_size_hmm,
        params.rfam_model_clans,
        params.rfam_models,
        ch_ssu_db,
        params.ssu_label,
        params.ssu_otus,
        params.ssu_tax,
        ch_egg_nog_data_dir,
        ch_egg_nog_db,
        ch_egg_nog_diamond_db
    )

    AFTER_QC(
        params.5.8s_pattern,
        params.5s_pattern,
        params.cgc_config,
        params.cgc_postfixes,
        params.hmm_gathering_bit_score,
        params.hmm_name_database,
        params.hmm_omit_alignment,
        params.inter_pro_scan_applications,
        params.inter_pro_scan_databases,
        params.inter_pro_scan_output_format,
        params.cgc_chunk_size,
        BEFORE_QC.out.filtered_fasta,
        params.func_ann_names_hmmer,
        params.func_ann_names_ips,
        params.go_config,
        params.hmmsearch_header,
        params.ips_header,
        params.ko_file,
        ch_lsu_db,
        params.lsu_label,
        params.lsu_otus,
        params.lsu_tax,
        BEFORE_QC.out.motus_input,
        params.other_ncrna_models,
        params.protein_chunk_size_ips,
        params.protein_chunk_size_hmm,
        params.rfam_model_clans,
        params.rfam_models,
        ch_ssu_db,
        params.ssu_label,
        params.ssu_otus,
        params.ssu_tax,
        ch_egg_nog_data_dir,
        ch_egg_nog_db,
        ch_egg_nog_diamond_db
    )

    BEFORE_QC(
        params.qc_min_length,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    BEFORE_QC(
        params.qc_min_length,
        ch_forward_reads,
        ch_reverse_reads,
        ch_single_reads
    )

    TOUCH_FILE_FLAG()

    TOUCH_FILE_FLAG()

    TOUCH_NO_CDS_FLAG()


}


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


ch_mediadb = Channel.fromPath( params.mediadb ).ifEmpty( null )
ch_bins    = Channel.fromPath( params.bins ).toList()


workflow  {

    CARVEME(
        PRODIGAL.out.predicted_proteins_faa,
        ch_mediadb
    )

    CARVEME_FILES_TO_FOLDER(
        COMPRESS_CARVEME.out.outfile.toList()
    )

    COMPRESS_CARVEME(
        CARVEME.out.carveme_gem
    )

    COMPRESS_PRODIGAL(
        PRODIGAL.out.predicted_proteins_faa
    )

    GEMSTATS(
        CARVEME.out.carveme_gem.toList()
    )

    MEMOTE_FILES_TO_FOLDER(
        MEMOTE_REPORT_SNAPSHOT.out.report_html.toList()
    )

    MEMOTE_REPORT_SNAPSHOT(
        CARVEME.out.carveme_gem
    )

    MEMOTE_RUN(
        CARVEME.out.carveme_gem
    )

    PRODIGAL(
        ch_bins.flatten().first()
    )

    PRODIGAL_FILES_TO_FOLDER(
        COMPRESS_PRODIGAL.out.outfile.toList()
    )

    SMETANA(
        CARVEME.out.carveme_gem.toList()
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


ch_metadata      = Channel.fromPath( params.metadata ).ifEmpty( null )
ch_reference_db  = Channel.fromPath( params.reference_db ).ifEmpty( null )
ch_reverse_reads = Channel.fromPath( params.reverse_reads ).ifEmpty( null )
ch_forward_reads = Channel.fromPath( params.forward_reads )


workflow  {

    FASTQC(
        ch_forward_reads.toList()
    )

    FASTQC_FILES_TO_FOLDER(
        FASTQC.out.html_files
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )

    READS_TO_FOLDER(
        ch_forward_reads.toList()
    )


}


ch_metadata      = Channel.fromPath( params.metadata ).ifEmpty( null )
ch_reference_db  = Channel.fromPath( params.reference_db ).ifEmpty( null )
ch_reverse_reads = Channel.fromPath( params.reverse_reads ).ifEmpty( null )
ch_forward_reads = Channel.fromPath( params.forward_reads )


workflow  {

    FASTQC(
        ch_forward_reads.toList()
    )

    FASTQC(
        ch_forward_reads.toList()
    )

    FASTQC_FILES_TO_FOLDER(
        FASTQC.out.html_files
    )

    FASTQC_FILES_TO_FOLDER(
        FASTQC.out.html_files
    )

    FOLDER_COMPRESSION(
        PICRUST2.out.intermediate
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX(
        ch_reference_db,
        READS_TO_FOLDER.out.results
    )

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_FILES_TO_FOLDER(
        NGTAX.out.biom.toList()
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    NGTAX_TO_TSV_FASTA(
        NGTAX.out.turtle,
        ch_metadata
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )

    PHYLOSEQ_FILES_TO_FOLDER(
        NGTAX_TO_TSV_FASTA.out.physeq_asv.toList()
    )

    PICRUST2(
        NGTAX_TO_TSV_FASTA.out.picrust_fasta,
        NGTAX_TO_TSV_FASTA.out.picrust_tsv
    )

    PICRUST_FILES_TO_FOLDER(
        PICRUST2.out.EC_predicted.tsv.gz.toList(),
        PICRUST2.out.EC_metagenome_out.toList()
    )

    READS_TO_FOLDER(
        ch_forward_reads.toList()
    )

    READS_TO_FOLDER(
        ch_forward_reads.toList()
    )


}


ch_assembly               = Channel.fromPath( params.assembly )
ch_illumina_forward_reads = Channel.fromPath( params.illumina_forward_reads )
ch_illumina_reverse_reads = Channel.fromPath( params.illumina_reverse_reads )


workflow  {

    BAM_INDEX(
        SAM_TO_SORTED_BAM.out.sortedbam
    )

    EXPRESSIONTOOL_BAM_INDEX(
        SAM_TO_SORTED_BAM.out.sortedbam,
        BAM_INDEX.out.bam_index
    )

    PILON(
        ch_assembly,
        EXPRESSIONTOOL_BAM_INDEX.out.hybrid_bamindex
    )

    READMAPPING_PILON(
        ch_illumina_forward_reads,
        ch_assembly,
        ch_illumina_reverse_reads
    )

    SAM_TO_SORTED_BAM(
        READMAPPING_PILON.out.sam
    )

    VCF_COMPRESS(
        PILON.out.pilon_vcf
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
