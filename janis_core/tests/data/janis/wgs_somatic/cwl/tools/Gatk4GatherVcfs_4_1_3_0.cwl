#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: Gather VCFs'
doc: |-
  GatherVcfs (Picard)
              
  Gathers multiple VCF files from a scatter operation into a single VCF file. 
  Input files must be supplied in genomic order and must not have events at overlapping positions.

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
- id: vcfs
  label: vcfs
  doc: '[default: []] (-I) Input VCF file(s).'
  type:
    type: array
    inputBinding:
      prefix: --INPUT
    items: File
  inputBinding: {}
- id: outputFilename
  label: outputFilename
  doc: '[default: null] (-O) Output VCF file.'
  type:
  - string
  - 'null'
  default: generated.gathered.vcf
  inputBinding:
    prefix: --OUTPUT
- id: argumentsFile
  label: argumentsFile
  doc: '[default: []] read one or more arguments files and add them to the command
    line'
  type:
  - type: array
    items: File
  - 'null'
  inputBinding:
    prefix: --arguments_file
- id: compressionLevel
  label: compressionLevel
  doc: |-
    [default: 5] Compression level for all compressed files created (e.g. BAM and VCF).
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --COMPRESSION_LEVEL
- id: createIndex
  label: createIndex
  doc: |-
    [default: TRUE] Whether to create a BAM index when writing a coordinate-sorted BAM file.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --CREATE_INDEX
- id: createMd5File
  label: createMd5File
  doc: |-
    [default: FALSE] Whether to create an MD5 digest for any BAM or FASTQ files created.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --CREATE_MD5_FILE
- id: ga4ghClientSecrets
  label: ga4ghClientSecrets
  doc: |-
    [default: client_secrets.json] Google Genomics API client_secrets.json file path.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --GA4GH_CLIENT_SECRETS
- id: maxRecordsInRam
  label: maxRecordsInRam
  doc: |-
    [default: 500000] When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --MAX_RECORDS_IN_RAM
- id: quiet
  label: quiet
  doc: '[default: FALSE] Whether to suppress job-summary info on System.err.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --QUIET
- id: referenceSequence
  label: referenceSequence
  doc: '[default: null] Reference sequence file.'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --REFERENCE_SEQUENCE
- id: tmpDir
  label: tmpDir
  doc: |-
    [default: []] One or more directories with space available to be used by this program for temporary storage of working files
  type: string
  default: /tmp
  inputBinding:
    prefix: --TMP_DIR
- id: useJdkDeflater
  label: useJdkDeflater
  doc: |-
    [default: FALSE] (-use_jdk_deflater) Use the JDK Deflater instead of the Intel Deflater for writing compressed output
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --USE_JDK_DEFLATER
- id: useJdkInflater
  label: useJdkInflater
  doc: |-
    [default: FALSE] (-use_jdk_inflater) Use the JDK Inflater instead of the Intel Inflater for reading compressed input
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --USE_JDK_INFLATER
- id: validationStringency
  label: validationStringency
  doc: |-
    [default: STRICT] Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --VALIDATION_STRINGENCY
- id: verbosity
  label: verbosity
  doc: '[default: INFO] Control verbosity of logging.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --VERBOSITY

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.gathered.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- GatherVcfs
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4GatherVcfs
