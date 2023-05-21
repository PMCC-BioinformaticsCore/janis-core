#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Fasta_to_Scaffolds2Bin

doc: Converts genome bins in fasta format to scaffolds-to-bin table. (DAS Tool helper script)

requirements:
  - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      dastool :
        version: ["1.1.4"]
        specs: ["https://anaconda.org/bioconda/das_tool"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/das_tool:1.1.4

baseCommand: [ "Fasta_to_Contig2Bin.sh" ]

inputs:
  binner_name:
    type: string
    doc: Binner name used to create the bins
    label: Binner name
  bin_folder:
    type: Directory
    doc: Input assembly in fasta format
    label: Input assembly
    inputBinding:
      prefix: --input_folder
  extension:
    type: string?
    doc: Extension of fasta files. (default fasta)
    label: Fasta extension
    inputBinding:
      prefix: --extension

stdout: $(inputs.binner_name)_Contig2Bin.tsv

outputs:
  table:
    type: File
    outputBinding:
      glob: $(inputs.binner_name)_Contig2Bin.tsv

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
s:dateCreated: "2022-09-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
