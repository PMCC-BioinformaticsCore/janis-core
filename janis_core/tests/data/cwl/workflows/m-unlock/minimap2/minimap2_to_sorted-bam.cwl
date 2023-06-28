#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

label: "Minimap2 to sorted bam"
doc: |
  Minimap2 directly piped to samtools for a sorted bam file
  minimap2 | samtools view | samtools sort

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

hints:
  SoftwareRequirement:
    packages:
      minimap2:
        version: ["2.24"]
        specs: ["https://anaconda.org/bioconda/minimap2"]
      samtools:
        version: ["1.15.1"]
        specs: ["https://anaconda.org/bioconda/samtools"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/minimap2:2.24

arguments:
    - shellQuote: false
      valueFrom: >
        ${
          var minimap2_path = "/unlock/infrastructure/binaries/minimap2/minimap2_v2.24/minimap2"
          var samtools_path = "/unlock/infrastructure/binaries/samtools/samtools_v1.15/bin/samtools"
          // minimap | samtools view | samtools sort
          var cmd = minimap2_path + " --split-prefix temp -a -t " + inputs.threads + " -x " + inputs.preset + " " + inputs.reference.path + " " + inputs.reads.path + "\
              | " + samtools_path + " view -@ " + inputs.threads + " -hu - \
              | " + samtools_path + " sort -@ " + inputs.threads + " -o " + inputs.identifier +".bam";
          return cmd;
        }

inputs:
  threads:
    label: Number of CPU-threads used by minimap2.
    type: int?
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  reference:
    label: Target sequence in FASTQ/FASTA format (can be gzipped).
    type: File
  reads:
    label: Query sequence in FASTQ/FASTA format (can be gzipped).
    type: File
  preset:
    label: read type
    doc: |
      - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
      - map-hifi - PacBio HiFi reads vs reference mapping
      - ava-pb/ava-ont - PacBio/Nanopore read overlap
      - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
      - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
      - sr - genomic short-read mapping
    type: string

outputs:
  bam: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).bam

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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
