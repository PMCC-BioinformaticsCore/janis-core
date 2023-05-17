#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "STAR Spliced Transcripts Alignment to a Reference"

doc: |
    Runs STAR in alignment mode

requirements:
 - class: InlineJavascriptRequirement
 - class: InitialWorkDirRequirement
   listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "STAR"
      writable: true

hints:
  SoftwareRequirement:
    packages:
      star:
        version: ["2.7.10a"]
        specs: ["https://anaconda.org/bioconda/star"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/star:2.7.10a

# baseCommand: [/unlock/infrastructure/binaries/STAR_v2.7.3a/bin/Linux_x86_64/STAR, --runMode, alignReads]   
baseCommand: [STAR, --runMode, alignReads]

inputs:
  genomeDir:
    type: Directory
    inputBinding:
      prefix: "--genomeDir"

  forward_reads:
    type:
     - File
     - File[]
    inputBinding:
      prefix: "--readFilesIn "
      separate: false
      itemSeparator: ","
      position: 1

  reverse_reads:
    type:
     - "null"
     - File
     - File[]
    inputBinding:
      prefix: ""
      separate: false
      itemSeparator: ","
      position: 2

  # Optional Inputs
  threads:
    type: int?
    inputBinding:
      prefix: "--runThreadN"
      
  OutFileNamePrefix:
    type: string?
    inputBinding:
      prefix: "--outFileNamePrefix"

  quantMode:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - TranscriptomeSAM
        - GeneCounts
    doc: Run with get gene quantification
    inputBinding:
      prefix: "--quantMode"

  sjdbGTFfile:
    type: File?
    inputBinding:
      prefix: "--sjdbGTFfile"

  Overhang:
    type: int?
    inputBinding:
      prefix: "--sjdbOverhang"

  sjdbGTFtagExonParentGene:
    type: string?
    doc: GTF attribute name for parent gene ID (default gene_id)
    inputBinding:
      prefix: "--sjdbGTFtagExonParentGene"

  sjdbGTFtagExonParentGeneName:
    type: string?
    doc: GTF attrbute name for parent gene name
    inputBinding:
      prefix: "--sjdbGTFtagExonParentGeneName"

  sjdbGTFtagExonParentGeneType:
    type: string?
    doc: GTF attrbute name for parent gene type
    inputBinding:
      prefix: "--sjdbGTFtagExonParentGeneType"

  OutFilterType:
    type:
     - "null"
     - type: enum
       symbols:
        - Normal
        - BySJout
    inputBinding:
      prefix: "--outFilterType"

  OutFilterIntronMotifs:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - RemoveNoncanonical
        - RemoveNoncanonicalUnannotated
    inputBinding:
      prefix: "--outFilterIntronMotifs"
  
  outSAMtype:
    type:
      type: array
      items: string
    default: [BAM, SortedByCoordinate]
    inputBinding:
      prefix: --outSAMtype
    doc: |
      strings: type of SAM/BAM output
      1st word:
      BAM  ... output BAM without sorting
      SAM  ... output SAM without sorting
      None ... no SAM/BAM output
      2nd, 3rd:
      Unsorted           ... standard unsorted
      SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.
  
  ReadFilesCommand:
    type: string?
    inputBinding:
      prefix: "--readFilesCommand"
    default: zcat

  AlignIntronMin:
    type: int?
    inputBinding:
      prefix: "--alignIntronMin"
  
  AlignIntronMax:
    type: int?
    inputBinding:
      prefix: "--alignIntronMax"
  
  AlignMatesGapMax:
    type: int?
    inputBinding:
      prefix: "--alignMatesGapMax"

  AlignSJoverhangMin:
    type: int?
    inputBinding:
      prefix: "--alignSJoverhangMin"
  
  AlignSJDBoverhangMin:
    type: int?
    inputBinding:
      prefix: "--alignSJDBoverhangMin"
  
  SeedSearchStartLmax:
    type: int?
    inputBinding:
      prefix: "--seedSearchStartLmax"

  ChimOutType:
    type:
     - "null"
     - type: enum
       symbols:
        - Junctions
        - SeparateSAMold
        - WithinBAM
        - "WithinBAM HardClip"
        - "WithinBAM SoftClip"

  ChimSegmentMin:
    type: int?
    inputBinding:
      prefix: "--chimSegmentMin"
    
  ChimJunctionOverhangMin:
    type: int?
    inputBinding:
      prefix: "--chimJunctionOverhangMin"

  OutFilterMultimapNmax:
    type: int?
    inputBinding:
      prefix: "--outFilterMultimapNmax"
  
  OutFilterMismatchNmax:
    type: int?
    inputBinding:
      prefix: "--outFilterMismatchNmax"

  OutFilterMismatchNoverLmax:
    type: double?
    inputBinding:
      prefix: "--outFilterMismatchNoverLmax"

  OutReadsUnmapped:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - Fastx
    inputBinding:
      prefix: "--outReadsUnmapped"
  
  OutSAMstrandField:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - intronMotif
    inputBinding:
      prefix: "--outSAMstrandField"
  
  OutSAMunmapped:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - Within
        - "Within KeepPairs"
    inputBinding:
      prefix: "--outSAMunmapped"
  
  OutSAMmapqUnique:
    type: int?
    inputBinding:
      prefix: "--outSAMmapqUnique"
  
  OutSamMode:
    type: 
     - "null"
     - type: enum
       symbols:
        - None
        - Full
        - NoQS
    inputBinding:
      prefix: "--outSAMmode"
  
  LimitOutSAMoneReadBytes:
    type: int?
    inputBinding:
      prefix: "--limitOutSAMoneReadBytes"
  
  GenomeLoad:
    type:
     - "null"
     - type: enum
       symbols:
        - LoadAndKeep
        - LoadAndRemove
        - LoadAndExit
        - Remove
        - NoSharedMemory
    inputBinding:
      prefix: "--genomeLoad"

outputs:
  sam:
    type:
     - "null"
     - File
     - File[]
    outputBinding:
      glob: "*.sam"

  bam:
    type:
     - "null"
     - File
     - File[]
    outputBinding:
      glob: "*.bam"

  unmapped_reads:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: "*Unmapped*"

  reads_per_gene:
    type:
      - "null"
      - File
    outputBinding:
      glob: "*ReadsPerGene.out.tab"

  splice_junctions:
    type:
      - "null"
      - File
    outputBinding:
      glob: "*SJ.out.tab"

  log_file:
    type: File
    outputBinding:
      glob: "*Log.out"

  final_log_file:
    type: File
    outputBinding:
      glob: "*Log.final.out"

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
 s: http://schema.org/