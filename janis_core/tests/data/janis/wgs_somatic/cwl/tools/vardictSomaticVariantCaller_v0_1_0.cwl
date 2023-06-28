#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.0
label: Vardict Somatic Variant Caller

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: normal_bam
  type: File
  secondaryFiles:
  - .bai
- id: tumor_bam
  type: File
  secondaryFiles:
  - .bai
- id: normal_name
  type: string
- id: tumor_name
  type: string
- id: intervals
  type: File
- id: allele_freq_threshold
  type: float
  default: 0.05
- id: header_lines
  type: File
- id: reference
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
- id: vardict_chromNamesAreNumbers
  doc: Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
  type: boolean
  default: true
- id: vardict_vcfFormat
  doc: VCF format output
  type: boolean
  default: true
- id: vardict_chromColumn
  doc: The column for chromosome
  type: int
  default: 1
- id: vardict_regStartCol
  doc: The column for region start, e.g. gene start
  type: int
  default: 2
- id: vardict_geneEndCol
  doc: The column for region end, e.g. gene end
  type: int
  default: 3
- id: compressvcf_stdout
  doc: 'c: Write to standard output, keep original files unchanged.'
  type: boolean
  default: true

outputs:
- id: variants
  type: File
  secondaryFiles:
  - .tbi
  outputSource: tabixvcf/out
- id: out
  type: File
  outputSource: filterpass/out

steps:
- id: vardict
  label: Vardict (Somatic)
  in:
  - id: tumorBam
    source: tumor_bam
  - id: normalBam
    source: normal_bam
  - id: intervals
    source: intervals
  - id: reference
    source: reference
  - id: tumorName
    source: tumor_name
  - id: normalName
    source: normal_name
  - id: alleleFreqThreshold
    source: allele_freq_threshold
  - id: chromNamesAreNumbers
    source: vardict_chromNamesAreNumbers
  - id: chromColumn
    source: vardict_chromColumn
  - id: geneEndCol
    source: vardict_geneEndCol
  - id: regStartCol
    source: vardict_regStartCol
  - id: vcfFormat
    source: vardict_vcfFormat
  run: vardict_somatic_1_6_0.cwl
  out:
  - id: out
- id: annotate
  label: 'BCFTools: Annotate'
  in:
  - id: vcf
    source: vardict/out
  - id: headerLines
    source: header_lines
  run: bcftoolsAnnotate_v1_5.cwl
  out:
  - id: out
- id: compressvcf
  label: BGZip
  in:
  - id: file
    source: annotate/out
  - id: stdout
    source: compressvcf_stdout
  run: bgzip_1_2_1.cwl
  out:
  - id: out
- id: tabixvcf
  label: Tabix
  in:
  - id: inp
    source: compressvcf/out
  run: tabix_1_2_1.cwl
  out:
  - id: out
- id: splitnormalisevcf
  label: Split Multiple Alleles
  in:
  - id: vcf
    source: annotate/out
  - id: reference
    source: reference
  run: SplitMultiAllele_v0_5772.cwl
  out:
  - id: out
- id: trim
  label: Trim IUPAC Bases
  in:
  - id: vcf
    source: splitnormalisevcf/out
  run: trimIUPAC_0_0_5.cwl
  out:
  - id: out
- id: filterpass
  label: Filter Vardict Somatic Vcf
  in:
  - id: vcf
    source: trim/out
  run: FilterVardictSomaticVcf_v1_9.cwl
  out:
  - id: out
id: vardictSomaticVariantCaller
