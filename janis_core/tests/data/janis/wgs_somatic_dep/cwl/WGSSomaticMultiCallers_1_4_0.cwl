#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: WGS Somatic (Multi callers)
doc: |
  This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs:

  - Takes raw sequence data in the FASTQ format;
  - align to the reference genome using BWA MEM;
  - Marks duplicates using Picard;
  - Call the appropriate somatic variant callers (GATK / Strelka / VarDict);
  - Outputs the final variants in the VCF format.

  **Resources**

  This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:

  - https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

  This pipeline expects the assembly references to be as they appear in that storage     (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
  The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: ScatterFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: MultipleInputFeatureRequirement

inputs:
- id: normal_inputs
  doc: |-
    An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
  type:
    type: array
    items:
      type: array
      items: File
- id: tumor_inputs
  doc: |-
    An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
  type:
    type: array
    items:
      type: array
      items: File
- id: normal_name
  doc: |-
    Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
  type: string
- id: tumor_name
  doc: |-
    Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
  type: string
- id: reference
  doc: |2-
        The reference genome from which to align the reads. This requires a number indexes (can be generated     with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

        This pipeline expects the assembly references to be as they appear in the GCP example. For example:
            - HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

        - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
- id: snps_dbsnp
  doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
  type: File
  secondaryFiles:
  - .tbi
- id: snps_1000gp
  doc: |-
    From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``. Accessible from the HG38 genomics-public-data google cloud bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/ 
  type: File
  secondaryFiles:
  - .tbi
- id: known_indels
  doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
  type: File
  secondaryFiles:
  - .tbi
- id: mills_indels
  doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
  type: File
  secondaryFiles:
  - .tbi
- id: gatk_intervals
  doc: List of intervals over which to split the GATK variant calling
  type:
    type: array
    items: File
- id: gridss_blacklist
  doc: |-
    BED file containing regions to ignore. For more information, visit: https://github.com/PapenfussLab/gridss#blacklist
  type: File
- id: vardict_intervals
  doc: List of intervals over which to split the VarDict variant calling
  type:
    type: array
    items: File
- id: strelka_intervals
  doc: An interval for which to restrict the analysis to.
  type: File
  secondaryFiles:
  - .tbi
- id: gnomad
  doc: |-
    The genome Aggregation Database (gnomAD). This VCF must be compressed and tabix indexed. This is specific for your genome (eg: hg38 / br37) and can usually be found with your reference. For example for HG38, the Broad institute provide the following af-only-gnomad compressed and tabix indexed VCF: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=af-only
  type: File
  secondaryFiles:
  - .tbi
- id: panel_of_normals
  doc: VCF file of sites observed in normal.
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
- id: cutadapt_adapters
  doc: |2-
                    Specifies a containment list for cutadapt, which contains a list of sequences to determine valid
                    overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets
                    of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
  type: File
- id: allele_freq_threshold
  doc: "The threshold for VarDict's allele frequency, default: 0.05 or 5%"
  type: float
  default: 0.05
- id: combine_variants_type
  doc: germline | somatic
  type: string
  default: somatic
- id: combine_variants_columns
  doc: Columns to keep, seperated by space output vcf (unsorted)
  type:
    type: array
    items: string
  default:
  - AD
  - DP
  - GT

outputs:
- id: out_normal_fastqc_reports
  type:
    type: array
    items:
      type: array
      items: File
  outputSource: normal/out_fastqc_reports
- id: out_tumor_fastqc_reports
  type:
    type: array
    items:
      type: array
      items: File
  outputSource: tumor/out_fastqc_reports
- id: out_normal_performance_summary
  doc: A text file of performance summary of NORMAL bam
  type: File
  outputSource: normal/out_performance_summary
- id: out_tumor_performance_summary
  doc: A text file of performance summary of TUMOR bam
  type: File
  outputSource: tumor/out_performance_summary
- id: out_normal_bam
  type: File
  secondaryFiles:
  - .bai
  outputSource: normal/out_bam
- id: out_tumor_bam
  type: File
  secondaryFiles:
  - .bai
  outputSource: tumor/out_bam
- id: out_gridss_assembly
  doc: Assembly returned by GRIDSS
  type: File
  outputSource: vc_gridss/assembly
- id: out_variants_gridss
  doc: Variants from the GRIDSS variant caller
  type: File
  outputSource: vc_gridss/out
- id: out_variants_gatk
  doc: Merged variants from the GATK caller
  type: File
  outputSource: vc_gatk_sort_combined/out
- id: out_variants_split
  doc: Unmerged variants from the GATK caller (by interval)
  type:
    type: array
    items: File
  outputSource: vc_gatk/out
- id: out_variants_vardict_split
  doc: Unmerged variants from the VarDict caller (by interval)
  type:
    type: array
    items: File
  outputSource: vc_vardict/out
- id: out_variants_vardict
  doc: Merged variants from the VarDict caller
  type: File
  outputSource: vc_vardict_sort_combined/out
- id: out_variants_strelka
  doc: Variants from the Strelka variant caller
  type: File
  outputSource: vc_strelka/out
- id: out_variants
  doc: Combined variants from GATK, VarDict and Strelka callers
  type: File
  outputSource: addbamstats/out

steps:
- id: tumor
  in:
  - id: reads
    source: tumor_inputs
  - id: sample_name
    source: tumor_name
  - id: reference
    source: reference
  - id: cutadapt_adapters
    source: cutadapt_adapters
  - id: gatk_intervals
    source: gatk_intervals
  - id: snps_dbsnp
    source: snps_dbsnp
  - id: snps_1000gp
    source: snps_1000gp
  - id: known_indels
    source: known_indels
  - id: mills_indels
    source: mills_indels
  run: tools/somatic_subpipeline.cwl
  out:
  - id: out_bam
  - id: out_fastqc_reports
  - id: out_performance_summary
- id: normal
  in:
  - id: reads
    source: normal_inputs
  - id: sample_name
    source: normal_name
  - id: reference
    source: reference
  - id: cutadapt_adapters
    source: cutadapt_adapters
  - id: gatk_intervals
    source: gatk_intervals
  - id: snps_dbsnp
    source: snps_dbsnp
  - id: snps_1000gp
    source: snps_1000gp
  - id: known_indels
    source: known_indels
  - id: mills_indels
    source: mills_indels
  run: tools/somatic_subpipeline.cwl
  out:
  - id: out_bam
  - id: out_fastqc_reports
  - id: out_performance_summary
- id: vc_gridss
  label: Gridss
  in:
  - id: bams
    source:
    - normal/out_bam
    - tumor/out_bam
  - id: reference
    source: reference
  - id: blacklist
    source: gridss_blacklist
  run: tools/gridss_v2_6_2.cwl
  out:
  - id: out
  - id: assembly
- id: bqsr_normal
  label: GATK Base Recalibration on Bam
  in:
  - id: bam
    source: normal/out_bam
  - id: intervals
    source: gatk_intervals
  - id: reference
    source: reference
  - id: snps_dbsnp
    source: snps_dbsnp
  - id: snps_1000gp
    source: snps_1000gp
  - id: known_indels
    source: known_indels
  - id: mills_indels
    source: mills_indels
  scatter:
  - intervals
  run: tools/GATKBaseRecalBQSRWorkflow_4_1_3.cwl
  out:
  - id: out
- id: bqsr_tumor
  label: GATK Base Recalibration on Bam
  in:
  - id: bam
    source: tumor/out_bam
  - id: intervals
    source: gatk_intervals
  - id: reference
    source: reference
  - id: snps_dbsnp
    source: snps_dbsnp
  - id: snps_1000gp
    source: snps_1000gp
  - id: known_indels
    source: known_indels
  - id: mills_indels
    source: mills_indels
  scatter:
  - intervals
  run: tools/GATKBaseRecalBQSRWorkflow_4_1_3.cwl
  out:
  - id: out
- id: vc_gatk
  label: GATK4 Somatic Variant Caller
  in:
  - id: normal_bam
    source: bqsr_normal/out
  - id: tumor_bam
    source: bqsr_tumor/out
  - id: normal_name
    source: normal_name
  - id: intervals
    source: gatk_intervals
  - id: reference
    source: reference
  - id: gnomad
    source: gnomad
  - id: panel_of_normals
    source: panel_of_normals
  scatter:
  - intervals
  - normal_bam
  - tumor_bam
  scatterMethod: dotproduct
  run: tools/GATK4_SomaticVariantCaller_4_1_3_0.cwl
  out:
  - id: variants
  - id: out_bam
  - id: out
- id: vc_gatk_merge
  label: 'GATK4: Gather VCFs'
  in:
  - id: vcfs
    source: vc_gatk/out
  run: tools/Gatk4GatherVcfs_4_1_3_0.cwl
  out:
  - id: out
- id: vc_gatk_compress_for_sort
  label: BGZip
  in:
  - id: file
    source: vc_gatk_merge/out
  run: tools/bgzip_1_2_1.cwl
  out:
  - id: out
- id: vc_gatk_sort_combined
  label: 'BCFTools: Sort'
  in:
  - id: vcf
    source: vc_gatk_compress_for_sort/out
  run: tools/bcftoolssort_v1_9.cwl
  out:
  - id: out
- id: vc_gatk_uncompress_for_combine
  label: UncompressArchive
  in:
  - id: file
    source: vc_gatk_sort_combined/out
  run: tools/UncompressArchive_v1_0_0.cwl
  out:
  - id: out
- id: addbamstats
  label: Annotate Bam Stats to Somatic Vcf Workflow
  in:
  - id: normal_id
    source: normal_name
  - id: tumor_id
    source: tumor_name
  - id: normal_bam
    source: normal/out_bam
  - id: tumor_bam
    source: tumor/out_bam
  - id: reference
    source: reference
  - id: vcf
    source: vc_gatk_uncompress_for_combine/out
  run: tools/AddBamStatsSomatic_v0_1_0.cwl
  out:
  - id: out
- id: generate_vardict_headerlines
  label: GenerateVardictHeaderLines
  in:
  - id: reference
    source: reference
  run: tools/GenerateVardictHeaderLines_v0_1_0.cwl
  out:
  - id: out
- id: vc_vardict
  label: Vardict Somatic Variant Caller
  in:
  - id: normal_bam
    source: normal/out_bam
  - id: tumor_bam
    source: tumor/out_bam
  - id: normal_name
    source: normal_name
  - id: tumor_name
    source: tumor_name
  - id: intervals
    source: vardict_intervals
  - id: allele_freq_threshold
    source: allele_freq_threshold
  - id: header_lines
    source: generate_vardict_headerlines/out
  - id: reference
    source: reference
  scatter:
  - intervals
  run: tools/vardictSomaticVariantCaller_v0_1_0.cwl
  out:
  - id: variants
  - id: out
- id: vc_vardict_merge
  label: 'GATK4: Gather VCFs'
  in:
  - id: vcfs
    source: vc_vardict/out
  run: tools/Gatk4GatherVcfs_4_1_3_0.cwl
  out:
  - id: out
- id: vc_vardict_compress_for_sort
  label: BGZip
  in:
  - id: file
    source: vc_vardict_merge/out
  run: tools/bgzip_1_2_1.cwl
  out:
  - id: out
- id: vc_vardict_sort_combined
  label: 'BCFTools: Sort'
  in:
  - id: vcf
    source: vc_vardict_compress_for_sort/out
  run: tools/bcftoolssort_v1_9.cwl
  out:
  - id: out
- id: vc_vardict_uncompress_for_combine
  label: UncompressArchive
  in:
  - id: file
    source: vc_vardict_sort_combined/out
  run: tools/UncompressArchive_v1_0_0.cwl
  out:
  - id: out
- id: vc_strelka
  label: Strelka Somatic Variant Caller
  in:
  - id: normal_bam
    source: normal/out_bam
  - id: tumor_bam
    source: tumor/out_bam
  - id: reference
    source: reference
  - id: intervals
    source: strelka_intervals
  run: tools/strelkaSomaticVariantCaller_v0_1_1.cwl
  out:
  - id: sv
  - id: variants
  - id: out
- id: combine_variants
  label: Combine Variants
  in:
  - id: vcfs
    source:
    - vc_gatk_uncompress_for_combine/out
    - vc_strelka/out
    - vc_vardict_uncompress_for_combine/out
  - id: type
    source: combine_variants_type
  - id: columns
    source: combine_variants_columns
  - id: normal
    source: normal_name
  - id: tumor
    source: tumor_name
  run: tools/combinevariants_0_0_8.cwl
  out:
  - id: out
- id: combined_compress
  label: BGZip
  in:
  - id: file
    source: combine_variants/out
  run: tools/bgzip_1_2_1.cwl
  out:
  - id: out
- id: combined_sort
  label: 'BCFTools: Sort'
  in:
  - id: vcf
    source: combined_compress/out
  run: tools/bcftoolssort_v1_9.cwl
  out:
  - id: out
- id: combined_uncompress
  label: UncompressArchive
  in:
  - id: file
    source: combined_sort/out
  run: tools/UncompressArchive_v1_0_0.cwl
  out:
  - id: out
- id: combined_addbamstats
  label: Annotate Bam Stats to Somatic Vcf Workflow
  in:
  - id: normal_id
    source: normal_name
  - id: tumor_id
    source: tumor_name
  - id: normal_bam
    source: normal/out_bam
  - id: tumor_bam
    source: tumor/out_bam
  - id: reference
    source: reference
  - id: vcf
    source: combined_uncompress/out
  run: tools/AddBamStatsSomatic_v0_1_0.cwl
  out:
  - id: out
id: WGSSomaticMultiCallers
