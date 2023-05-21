#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "MultiQC reporting"

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "medaka_output"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          source /root/miniconda/bin/activate
          conda init bash
          conda activate /unlock/infrastructure/conda/multiqc/multiqc_v1.12.0
          multiqc --filename multiqc_report $@
doc: |
    Generates a MultiQC report of a given folder

hints:
  SoftwareRequirement:
    packages:
      multiqc:
        version: ["1.13a"]
        specs: ["https://anaconda.org/bioconda/multiqc"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/multiqc:1.13a

# stdout: multiqc.log
# stderr: multiqc.error

inputs:
  folder:
      type: string
      doc: Directory path to be processed by MultiQC
      label: MultiQC directory
      inputBinding:
        position: 100

baseCommand: ["bash", "script.sh"]

# arguments: ["multiqc","--filename","multiqc_report"]

outputs:
  # info:
  #   type: stdout
  # error: 
  #   type: stderr
  multiqc: 
    type: File
    outputBinding:
      glob: multiqc_report.*

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
s:dateCreated: "2022-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
