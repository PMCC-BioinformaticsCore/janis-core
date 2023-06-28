version development

import "tools/somatic_subpipeline.wdl" as S
import "tools/gridss_v2_6_2.wdl" as G
import "tools/GATKBaseRecalBQSRWorkflow_4_1_3.wdl" as G2
import "tools/GATK4_SomaticVariantCaller_4_1_3_0.wdl" as G3
import "tools/Gatk4GatherVcfs_4_1_3_0.wdl" as G4
import "tools/bgzip_1_2_1.wdl" as B
import "tools/bcftoolssort_v1_9.wdl" as B2
import "tools/UncompressArchive_v1_0_0.wdl" as U
import "tools/AddBamStatsSomatic_v0_1_0.wdl" as A
import "tools/GenerateVardictHeaderLines_v0_1_0.wdl" as G5
import "tools/vardictSomaticVariantCaller_v0_1_0.wdl" as V
import "tools/strelkaSomaticVariantCaller_v0_1_1.wdl" as S2
import "tools/combinevariants_0_0_8.wdl" as C

workflow WGSSomaticMultiCallers {
  input {
    Array[Array[File]] normal_inputs
    Array[Array[File]] tumor_inputs
    String normal_name
    String tumor_name
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    File snps_dbsnp
    File snps_dbsnp_tbi
    File snps_1000gp
    File snps_1000gp_tbi
    File known_indels
    File known_indels_tbi
    File mills_indels
    File mills_indels_tbi
    Array[File] gatk_intervals
    File gridss_blacklist
    Array[File] vardict_intervals
    File strelka_intervals
    File strelka_intervals_tbi
    File gnomad
    File gnomad_tbi
    File? panel_of_normals
    File? panel_of_normals_tbi
    File cutadapt_adapters
    Float? allele_freq_threshold = 0.05
    String? combine_variants_type = "somatic"
    Array[String]? combine_variants_columns = ["AD", "DP", "GT"]
  }
  call S.somatic_subpipeline as tumor {
    input:
      reads=tumor_inputs,
      sample_name=tumor_name,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      cutadapt_adapters=cutadapt_adapters,
      gatk_intervals=gatk_intervals,
      snps_dbsnp=snps_dbsnp,
      snps_dbsnp_tbi=snps_dbsnp_tbi,
      snps_1000gp=snps_1000gp,
      snps_1000gp_tbi=snps_1000gp_tbi,
      known_indels=known_indels,
      known_indels_tbi=known_indels_tbi,
      mills_indels=mills_indels,
      mills_indels_tbi=mills_indels_tbi
  }
  call S.somatic_subpipeline as normal {
    input:
      reads=normal_inputs,
      sample_name=normal_name,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      cutadapt_adapters=cutadapt_adapters,
      gatk_intervals=gatk_intervals,
      snps_dbsnp=snps_dbsnp,
      snps_dbsnp_tbi=snps_dbsnp_tbi,
      snps_1000gp=snps_1000gp,
      snps_1000gp_tbi=snps_1000gp_tbi,
      known_indels=known_indels,
      known_indels_tbi=known_indels_tbi,
      mills_indels=mills_indels,
      mills_indels_tbi=mills_indels_tbi
  }
  call G.gridss as vc_gridss {
    input:
      bams=[normal.out_bam, tumor.out_bam],
      bams_bai=[normal.out_bam_bai, tumor.out_bam_bai],
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      blacklist=gridss_blacklist
  }
  scatter (g in gatk_intervals) {
     call G2.GATKBaseRecalBQSRWorkflow as bqsr_normal {
      input:
        bam=normal.out_bam,
        bam_bai=normal.out_bam_bai,
        intervals=g,
        reference=reference,
        reference_fai=reference_fai,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_sa=reference_sa,
        reference_dict=reference_dict,
        snps_dbsnp=snps_dbsnp,
        snps_dbsnp_tbi=snps_dbsnp_tbi,
        snps_1000gp=snps_1000gp,
        snps_1000gp_tbi=snps_1000gp_tbi,
        known_indels=known_indels,
        known_indels_tbi=known_indels_tbi,
        mills_indels=mills_indels,
        mills_indels_tbi=mills_indels_tbi
    }
  }
  scatter (g in gatk_intervals) {
     call G2.GATKBaseRecalBQSRWorkflow as bqsr_tumor {
      input:
        bam=tumor.out_bam,
        bam_bai=tumor.out_bam_bai,
        intervals=g,
        reference=reference,
        reference_fai=reference_fai,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_sa=reference_sa,
        reference_dict=reference_dict,
        snps_dbsnp=snps_dbsnp,
        snps_dbsnp_tbi=snps_dbsnp_tbi,
        snps_1000gp=snps_1000gp,
        snps_1000gp_tbi=snps_1000gp_tbi,
        known_indels=known_indels,
        known_indels_tbi=known_indels_tbi,
        mills_indels=mills_indels,
        mills_indels_tbi=mills_indels_tbi
    }
  }
  scatter (Q in zip(gatk_intervals, zip(transpose([bqsr_normal.out, bqsr_normal.out_bai]), transpose([bqsr_tumor.out, bqsr_tumor.out_bai])))) {
     call G3.GATK4_SomaticVariantCaller as vc_gatk {
      input:
        normal_bam=Q.right.left[0],
        normal_bam_bai=Q.right.left[1],
        tumor_bam=Q.right.right[0],
        tumor_bam_bai=Q.right.right[1],
        normal_name=normal_name,
        intervals=Q.left,
        reference=reference,
        reference_fai=reference_fai,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_sa=reference_sa,
        reference_dict=reference_dict,
        gnomad=gnomad,
        gnomad_tbi=gnomad_tbi,
        panel_of_normals=panel_of_normals,
        panel_of_normals_tbi=panel_of_normals_tbi
    }
  }
  call G4.Gatk4GatherVcfs as vc_gatk_merge {
    input:
      vcfs=vc_gatk.out
  }
  call B.bgzip as vc_gatk_compress_for_sort {
    input:
      file=vc_gatk_merge.out
  }
  call B2.bcftoolssort as vc_gatk_sort_combined {
    input:
      vcf=vc_gatk_compress_for_sort.out
  }
  call U.UncompressArchive as vc_gatk_uncompress_for_combine {
    input:
      file=vc_gatk_sort_combined.out
  }
  call A.AddBamStatsSomatic as addbamstats {
    input:
      normal_id=normal_name,
      tumor_id=tumor_name,
      normal_bam=normal.out_bam,
      normal_bam_bai=normal.out_bam_bai,
      tumor_bam=tumor.out_bam,
      tumor_bam_bai=tumor.out_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      vcf=vc_gatk_uncompress_for_combine.out
  }
  call G5.GenerateVardictHeaderLines as generate_vardict_headerlines {
    input:
      reference=reference,
      reference_dict=reference_dict
  }
  scatter (v in vardict_intervals) {
     call V.vardictSomaticVariantCaller as vc_vardict {
      input:
        normal_bam=normal.out_bam,
        normal_bam_bai=normal.out_bam_bai,
        tumor_bam=tumor.out_bam,
        tumor_bam_bai=tumor.out_bam_bai,
        normal_name=normal_name,
        tumor_name=tumor_name,
        intervals=v,
        allele_freq_threshold=select_first([allele_freq_threshold, 0.05]),
        header_lines=generate_vardict_headerlines.out,
        reference=reference,
        reference_fai=reference_fai,
        reference_amb=reference_amb,
        reference_ann=reference_ann,
        reference_bwt=reference_bwt,
        reference_pac=reference_pac,
        reference_sa=reference_sa,
        reference_dict=reference_dict
    }
  }
  call G4.Gatk4GatherVcfs as vc_vardict_merge {
    input:
      vcfs=vc_vardict.out
  }
  call B.bgzip as vc_vardict_compress_for_sort {
    input:
      file=vc_vardict_merge.out
  }
  call B2.bcftoolssort as vc_vardict_sort_combined {
    input:
      vcf=vc_vardict_compress_for_sort.out
  }
  call U.UncompressArchive as vc_vardict_uncompress_for_combine {
    input:
      file=vc_vardict_sort_combined.out
  }
  call S2.strelkaSomaticVariantCaller as vc_strelka {
    input:
      normal_bam=normal.out_bam,
      normal_bam_bai=normal.out_bam_bai,
      tumor_bam=tumor.out_bam,
      tumor_bam_bai=tumor.out_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      intervals=strelka_intervals,
      intervals_tbi=strelka_intervals_tbi
  }
  call C.combinevariants as combine_variants {
    input:
      vcfs=[vc_gatk_uncompress_for_combine.out, vc_strelka.out, vc_vardict_uncompress_for_combine.out],
      type=select_first([combine_variants_type, "somatic"]),
      columns=select_first([combine_variants_columns, ["AD", "DP", "GT"]]),
      normal=normal_name,
      tumor=tumor_name
  }
  call B.bgzip as combined_compress {
    input:
      file=combine_variants.out
  }
  call B2.bcftoolssort as combined_sort {
    input:
      vcf=combined_compress.out
  }
  call U.UncompressArchive as combined_uncompress {
    input:
      file=combined_sort.out
  }
  call A.AddBamStatsSomatic as combined_addbamstats {
    input:
      normal_id=normal_name,
      tumor_id=tumor_name,
      normal_bam=normal.out_bam,
      normal_bam_bai=normal.out_bam_bai,
      tumor_bam=tumor.out_bam,
      tumor_bam_bai=tumor.out_bam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_amb=reference_amb,
      reference_ann=reference_ann,
      reference_bwt=reference_bwt,
      reference_pac=reference_pac,
      reference_sa=reference_sa,
      reference_dict=reference_dict,
      vcf=combined_uncompress.out
  }
  output {
    Array[Array[File]] out_normal_fastqc_reports = normal.out_fastqc_reports
    Array[Array[File]] out_tumor_fastqc_reports = tumor.out_fastqc_reports
    File out_normal_performance_summary = normal.out_performance_summary
    File out_tumor_performance_summary = tumor.out_performance_summary
    File out_normal_bam = normal.out_bam
    File out_normal_bam_bai = normal.out_bam_bai
    File out_tumor_bam = tumor.out_bam
    File out_tumor_bam_bai = tumor.out_bam_bai
    File out_gridss_assembly = vc_gridss.assembly
    File out_variants_gridss = vc_gridss.out
    File out_variants_gatk = vc_gatk_sort_combined.out
    Array[File] out_variants_split = vc_gatk.out
    Array[File] out_variants_vardict_split = vc_vardict.out
    File out_variants_vardict = vc_vardict_sort_combined.out
    File out_variants_strelka = vc_strelka.out
    File out_variants = addbamstats.out
  }
}