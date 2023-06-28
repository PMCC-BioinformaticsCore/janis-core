#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

label: "Minimap2 to (un)mapped long reads"
doc: |
  Get unmapped or mapped long reads reads in fastq.gz format using minimap2 and samtools. Mainly used for contamination removal.
   - requires pigz!
  minimap2 | samtools | pigz

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "minimap_run"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          minimap2 -a -t $5 -x $4 $2 $3 | samtools fastq -@ $5 -n $1 4 | pigz -p $5 > $6_filtered.fastq.gz

baseCommand: [ bash, -x, script.sh ]

hints:
  SoftwareRequirement:
    packages:
      minimap2:
        version: ["2.24"]
        specs: ["https://anaconda.org/bioconda/minimap2"]
      samtools:
        version: ["1.15.1"]
        specs: ["https://anaconda.org/bioconda/samtools"]
      pigz:
        version: ["2.6"]
        specs: ["https://anaconda.org/conda-forge/pigz"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/minimap2:2.24

arguments:
  - |
    ${
      if (inputs.output_mapped){
        return '-F';
      } else {
        return '-f';
      }
    }


inputs:
  threads:
    type: int?
    doc: Number of CPU-threads used by minimap2.
    label: Threads
    default: 4
    inputBinding:
      position: 4
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier
    inputBinding:
      position: 5
  reference:
    type: File
    doc: Target sequence in FASTQ/FASTA format (can be gzipped).
    label: Reference
    inputBinding:
      position: 1

  reads:
    type: File
    doc: Query sequence in FASTQ/FASTA format (can be gzipped).
    label: Reads
    inputBinding:
      position: 2
    
  output_mapped:
    type: boolean?
    doc: Keep reads mapped to the reference (default = output unmapped)
    label: Keep mapped
    default: false
    inputBinding:
      position: 6

  preset:
    type: string
    doc: |
      - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
      - map-hifi - PacBio HiFi reads vs reference mapping
      - ava-pb/ava-ont - PacBio/Nanopore read overlap
      - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
      - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
      - sr - genomic short-read mapping
    label: Read type
    inputBinding:
      position: 3

stderr: $(inputs.identifier)_minimap2.log

outputs:
  fastq: 
    type: File
    outputBinding:
      glob: $(inputs.identifier)_filtered.fastq.gz
  log: 
    type: stderr

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-5516-8391
    s:email: mailto:german.royvalgarcia@wur.nl
    s:name: Germán Royval
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
s:dateCreated: "2022-03-00"
s:dateModified: "2022-04-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/