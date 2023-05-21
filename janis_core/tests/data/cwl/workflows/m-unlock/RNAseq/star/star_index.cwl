#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "STAR indexer"

doc: |
    Creates the Genome index for STAR spliced RNAseq aligner (single fasta)

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

inputs:
  inputFile:
    type: File
    inputBinding:
      prefix: "--genomeFastaFiles"

#  indexFolder:
#    type: string
#    inputBinding:
#      prefix: "--genomeDir"

# Optional Inputs
  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
 
  sjdbGTFfile:
    type: File?
    inputBinding:
      prefix: "--sjdbGTFfile"

  sjdbOverhang:
    type: int?
    inputBinding:
      prefix: "--sjdbOverhang"
  
  sjdbFileChrStartEnd:
    type: string?
    inputBinding:
      prefix: "--sjdbFileChrStartEnd"
  
  genomeSAindexNbases:
    type: int?
    inputBinding:
      prefix: "--genomeSAindexNbases"

  genomeChrBinNbits:
    type: int?
    inputBinding:
      prefix: "--genomeChrBinNbits"

# baseCommand: [/unlock/infrastructure/binaries/STAR-2.7.3a/bin/Linux_x86_64_static/STAR, "--runMode", genomeGenerate]
baseCommand: [STAR, --runMode, genomeGenerate]

arguments: ["--genomeDir", STAR]

outputs:
  STAR: 
    type: Directory
    outputBinding:
      glob: "STAR"


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
  s: https://schema.org/
