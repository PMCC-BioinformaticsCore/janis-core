#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "aggregateBinDepths"

doc: |
    Aggregate bin depths using MetaBat2 using the script aggregateBinDepths.pl

requirements:
 - class: InlineJavascriptRequirement
hints:
  SoftwareRequirement:
    packages:
      metabat2:
        version: ["2.15"]
        specs: ["https://anaconda.org/bioconda/metabat2"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/metabat2:2.15

inputs:
  identifier:
    type: string
    doc: Name of the output file
    label: output file name

  metabatdepthsFile:
    type: File
    doc: Contig depths files obtained from metabat2 script jgi_summarize_bam_contig_depths
    label: contigs depths
    inputBinding:
      position: 1
  bins:
    type: File[]
    doc: Bin fasta files
    label: Bin fasta files
    inputBinding:
      position: 2

# baseCommand: [/unlock/infrastructure/binaries/MetaBAT/metabat_v2.12.1/aggregateBinDepths.pl]
baseCommand: [aggregateBinDepths.pl]

stdout: $(inputs.identifier)_binDepths.tsv

outputs:
  binDepths:
    type: File
    outputBinding:
        glob: $(inputs.identifier)_binDepths.tsv

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
