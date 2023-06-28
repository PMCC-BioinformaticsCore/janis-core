#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: HUMAnN Analysis

doc: |
  Runs the HUMAnN meta-omics taxonomic and functional profiling tool.

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      humann2:
        version: ["2.8.1"]
        specs: ["https://anaconda.org/bioconda/humann2"]      
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/humann2:2.8.1

inputs:
  fasta:
    type: File
    doc: FASTA of unaligned sequences
    label: Input fasta
    inputBinding:
      prefix: -i

  metaphlan_db:
    type: string
    doc: location of a indexed metaphlan database
    label: metaphlan database
    default: "--bowtie2db /unlock/references/databases/HUMAnN/metaphlan_databases/"
    inputBinding:
      prefix: --metaphlan-options

  threads:
    type: int?
    default: 2
    inputBinding:
      prefix: --threads
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

baseCommand: [humann]

arguments:
  - prefix: "-o"
    valueFrom: $(inputs.identifier)_HUMAnN
  - prefix: "--o-log"
    valueFrom: $(inputs.identifier)_HUMAnN.log

stdout: $(inputs.identifier)_HUMAnN.stdout.log

outputs:
  genefamilies_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_HUMAnN/$(inputs.identifier)_genefamilies.tsv
  pathabundance_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_HUMAnN/$(inputs.identifier)_pathabundance.tsv
  pathcoverage_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_HUMAnN/$(inputs.identifier)_pathcoverage.tsv
 
  log_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_HUMAnN.log
  
  stdout_out:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_HUMAnN.stdout.log

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
