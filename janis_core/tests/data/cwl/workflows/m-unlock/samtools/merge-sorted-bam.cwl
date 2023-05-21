#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "sam to sorted bam"
doc: |
  samtools merge out in

requirements:
 - class: ShellCommandRequirement
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      samtools:
        version: ["1.15.1"]
        specs: ["https://anaconda.org/bioconda/samtools"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  bam:
    type: File[]
    doc: sorted bam files
    label: sorted bam files

  threads:
      type: int
      doc: number of cpu threads to use
      label: cpu threads
      default: 1

outputs:
  sortedbam: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).sorted.bam

# samtools merge out in
arguments:
    - shellQuote: false
      valueFrom: >
        ${
          // var samtools_path = "/unlock/infrastructure/binaries/samtools/samtools_v1.15/bin/samtools"
          var samtools_path = "samtools"
          var cmd = samtools_path + " merge "+inputs.identifier+".sorted.bam "+ inputs.bam.map( function(x) {return  x.path} ).join(' ');
          return cmd;
        }

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
s:dateModified: "2022-02-22"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
