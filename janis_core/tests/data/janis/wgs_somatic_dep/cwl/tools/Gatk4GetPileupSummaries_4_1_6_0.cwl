#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: GetPileupSummaries'
doc: |-
  Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
  The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.6.0

inputs:
- id: javaOptions
  label: javaOptions
  type:
  - type: array
    items: string
  - 'null'
- id: compression_level
  label: compression_level
  doc: |-
    Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
  type:
  - int
  - 'null'
- id: bam
  label: bam
  doc: The SAM/BAM/CRAM file containing reads.
  type:
    type: array
    inputBinding:
      prefix: -I
    items: File
  inputBinding:
    position: 0
- id: sites
  label: sites
  doc: sites of common biallelic variants
  type: File
  secondaryFiles:
  - .tbi
  inputBinding:
    prefix: -V
- id: intervals
  label: intervals
  doc: -L (BASE) One or more genomic intervals over which to operate
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --intervals
- id: pileupTableOut
  label: pileupTableOut
  type:
  - string
  - 'null'
  default: generated.txt
  inputBinding:
    prefix: -O
    position: 1
- id: reference
  label: reference
  doc: reference to use when decoding CRAMS
  type:
  - File
  - 'null'
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: -R

outputs:
- id: out
  label: out
  doc: Table containing the pileup info
  type: File
  outputBinding:
    glob: generated.txt
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- GetPileupSummaries
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 64, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4GetPileupSummaries
