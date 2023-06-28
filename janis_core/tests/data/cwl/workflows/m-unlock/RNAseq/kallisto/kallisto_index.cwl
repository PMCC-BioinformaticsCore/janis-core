#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "kallisto indexer"

doc: |-
    Creates the index for pseudoalignment tool kallisto
    https://github.com/common-workflow-library/bio-cwl-tools/tree/release/Kallisto

requirements:
 - class: InlineJavascriptRequirement
 - class: InitialWorkDirRequirement
   listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "kallisto"
      writable: true

hints:
  SoftwareRequirement:
    packages:
      kallisto:
        version: ["0.48.0"]
        specs: ["https://anaconda.org/bioconda/kallisto"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/kallisto:0.48.0

inputs:
  inputFile:
    type: File
    inputBinding:
      position: 200
    
#Optional arguments
  kmerSize:
    type: int?
    inputBinding:
      prefix: "--kmer-size="
      separate: false

  makeUnique:
    type: boolean?
    inputBinding:
      prefix: "--make-unique"

arguments:
  - prefix: "--index="
    separate: false
    valueFrom: $("kallisto/" + inputs.inputFile.nameroot)_kallisto.idx

# baseCommand: [/unlock/infrastructure/binaries/kallisto/kallisto_v0.48.0/kallisto, index]
baseCommand: [kallisto, index]

outputs:
  kallisto_indexFolder: 
    type: Directory
    outputBinding:
      glob: kallisto

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
