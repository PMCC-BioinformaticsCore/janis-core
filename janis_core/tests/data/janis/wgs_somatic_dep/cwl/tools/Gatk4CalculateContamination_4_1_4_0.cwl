#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: CalculateContamination'
doc: |-
  Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.

  This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 for a step-by-step description of the workflow and Article#11127 for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.

  This tool borrows from ContEst by Cibulskis et al the idea of estimating contamination from ref reads at hom alt sites. However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number variation and independent contaminating reads. That is, ContEst assumes that each contaminating read is drawn randomly and independently from a different human. This tool uses a simpler estimate of contamination that relaxes these assumptions. In particular, it works in the presence of copy number variations and with an arbitrary number of contaminating samples. In addition, this tool is designed to work well with no matched normal data. However, one can run GetPileupSummaries on a matched normal bam file and input the result to this tool.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.4.0

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
- id: contaminationTable
  label: contaminationTable
  doc: Tables containing contamination information.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --contamination-table
- id: statsFile
  label: statsFile
  doc: The Mutect stats file output by Mutect2
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --stats
- id: readOrientationModel
  label: readOrientationModel
  doc: |-
    One or more .tar.gz files containing tables of prior artifact probabilities for the read orientation filter model, one table per tumor sample
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --orientation-bias-artifact-priors
- id: pileupTable
  label: pileupTable
  doc: pileup table from summarize pileup
  type: File
  inputBinding:
    prefix: -I
- id: segmentationFileOut
  label: segmentationFileOut
  doc: Reference sequence file
  type:
  - string
  - 'null'
  default: generated.mutect2_segments
  inputBinding:
    prefix: --tumor-segmentation
    valueFrom: $(inputs.pileupTable.basename).mutect2_segments
- id: contaminationFileOut
  label: contaminationFileOut
  type:
  - string
  - 'null'
  default: generated.mutect2_contamination
  inputBinding:
    prefix: -O
    position: 2
    valueFrom: $(inputs.pileupTable.basename).mutect2_contamination

outputs:
- id: contOut
  label: contOut
  doc: contamination Table
  type: File
  outputBinding:
    glob: $(inputs.pileupTable.basename).mutect2_contamination
    loadContents: false
- id: segOut
  label: segOut
  doc: segmentation based on baf
  type: File
  outputBinding:
    glob: $(inputs.pileupTable.basename).mutect2_segments
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- CalculateContamination
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4CalculateContamination
