#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Unbinned contigs"

doc: |
    Get unbinned contigs of the assembbly from a set of binned fasta files in fasta format (compressed).

requirements:
 - class: InlineJavascriptRequirement 

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  threads:
    type: int?
    label: Number of threads to use
    default: 8
  
  assembly_fasta:
    type: File
    doc: fasta file that was used for the binning.
    label: assembly fasta file

  bin_dir:
    type: Directory
    doc: folder containing bins in fasta format from metagenomic binning
    label: bins folder

arguments:
  - "$(inputs.assembly_fasta)"
  - "$(inputs.bin_dir)"
  - "$(inputs.identifier)"
  - "$(inputs.threads)"

baseCommand: [bash, "/unlock/infrastructure/scripts/metagenomics/get_unbinned_contigs.sh"]

outputs:
  unbinned_fasta:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_unbinned.fasta.gz

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
s:dateCreated: "2021-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
