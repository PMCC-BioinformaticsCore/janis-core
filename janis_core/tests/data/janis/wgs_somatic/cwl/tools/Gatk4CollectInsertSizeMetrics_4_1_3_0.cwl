#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: CollectInsertSizeMetrics'
doc: |-
  Provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.3.0

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
  doc: Input SAM or BAM file.  Required.
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: -I
    position: 10
- id: outputFilename
  label: outputFilename
  doc: File to write the output to.  Required.
  type:
  - string
  - 'null'
  default: generated.metrics.txt
  inputBinding:
    prefix: -O
- id: outputHistogram
  label: outputHistogram
  doc: 'File to write insert size Histogram chart to.  Required. '
  type:
  - string
  - 'null'
  default: generated.histogram.pdf
  inputBinding:
    prefix: -H
- id: argumentsFile
  label: argumentsFile
  doc: read one or more arguments files and add them to the command line
  type:
  - type: array
    inputBinding:
      prefix: --arguments_file
    items: File
  - 'null'
  inputBinding:
    position: 10
- id: assumeSorted
  label: assumeSorted
  doc: |-
    If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --ASSUME_SORTED
    position: 11
- id: deviations
  label: deviations
  doc: |-
    Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0. 
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --DEVIATIONS
    position: 11
- id: histogramWidth
  label: histogramWidth
  doc: |-
    Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --HISTOGRAM_WIDTH
    position: 11
- id: includeDuplicates
  label: includeDuplicates
  doc: |-
    If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --INCLUDE_DUPLICATES
    position: 11
- id: metricAccumulationLevel
  label: metricAccumulationLevel
  doc: |-
    The level(s) at  which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} .
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --METRIC_ACCUMULATION_LEVEL
    position: 11
- id: minimumPCT
  label: minimumPCT
  doc: |-
    When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05.
  type:
  - float
  - 'null'
  inputBinding:
    prefix: --MINIMUM_PCT
    position: 11
- id: stopAfter
  label: stopAfter
  doc: 'Stop after  processing N reads, mainly for debugging.  Default value: 0. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --STOP_AFTER
    position: 11
- id: version
  label: version
  doc: |-
    display the version number for this tool Default value: false. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --version
    position: 11
- id: showHidden
  label: showHidden
  doc: |-
    display hidden  arguments  Default  value: false.  Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --showHidden
    position: 11

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.metrics.txt
    loadContents: false
- id: outHistogram
  label: outHistogram
  type: File
  outputBinding:
    glob: generated.histogram.pdf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- CollectInsertSizeMetrics
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4CollectInsertSizeMetrics
