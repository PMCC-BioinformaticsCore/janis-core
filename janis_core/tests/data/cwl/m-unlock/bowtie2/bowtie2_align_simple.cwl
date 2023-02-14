#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Bowtie2 alignment"

doc: |
    Align reads to indexed genome. Stripped simple version; only paired end reads and sam output.

requirements:
 - class: InlineJavascriptRequirement
 - class: InitialWorkDirRequirement
   listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "bowtie2"
      writable: true
hints:
  SoftwareRequirement:
    packages:
      bowtie2:
        version: ["2.4.5"]
        specs: ["https://anaconda.org/bioconda/bowtie2"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/bowtie2:2.4.5

inputs:
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: --threads

  identifier:
    type: string
    doc: prefix of the filename outputs      
    default: "out"

  indexfolder:
    type: Directory?
    doc: "Folder with indices files"
    inputBinding:
      prefix: '-x'
      valueFrom: |
        ${
            for (var i = 0; i < self.listing.length; i++) {
                if (self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.bt2' ||
                    self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.bt2l'){
                  return self.listing[i].path.split('.').slice(0,-3).join('.');
                }
            }
            return null;
        } 

  forward_reads:
    type: File
    inputBinding:
      prefix: "-1"

  reverse_reads:
    type: File
    inputBinding:
      prefix: "-2"

arguments:
  - prefix: "-S"
    valueFrom: $(inputs.identifier)_bowtie2.sam
  - prefix: "--met-file"
    valueFrom: $(inputs.identifier)_bowtie2_metrics.txt

baseCommand: [/unlock/infrastructure/binaries/bowtie2/bowtie2_v2.4.5/bowtie2]
# baseCommand: [bowtie2]

outputs:
  sam: 
    type: File
    outputBinding:
      glob: "*.sam"
  metricsfile:
    type: File
    outputBinding:
      glob: "*metrics.txt"

  bowtie2_stats:
    type: stderr
stderr: $(inputs.identifier)_bowtie2_stats.txt

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
s:dateModified: "2022-02-23"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
