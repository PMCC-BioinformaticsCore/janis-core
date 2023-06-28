#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Gene prediction"

doc: |
    Runs KEGG KO annotation on protein sequences using SAPP. Requires a kofamscan conda environment.

hints:
  SoftwareRequirement:
    packages:
      java:
        version: ["17.0.3"]
        specs: ["https://anaconda.org/conda-forge/openjdk"]
      hmmer:
        version: ["3.2.1"]
        specs: ["https://anaconda.org/bioconda/hmmer"]
      kofamscan:
        version: ["1.3.0"]
        specs: ["https://anaconda.org/bioconda/kofamscan"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/kofamscan:1.3.0

requirements:
  InlineJavascriptRequirement: {}
  # InitialWorkDirRequirement:
    # listing:
      # - entryname: script.sh
        # entry: |-
          #!/bin/bash
          # source /root/miniconda/bin/activate
          # conda init bash
          # conda activate /unlock/infrastructure/conda/kofamscan_v1.3.0
          # java -Xmx5g -jar /SAPP-2.0.jar -kofamscan $@

inputs:
  input:
    type: File
    doc: Reference genome file used in RDF format
    label: Reference genome
    inputBinding:
      prefix: -input
  identifier:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
  threads:
    type: int?
    default: 3
    inputBinding:
      prefix: -threads

baseCommand: ["java", "-Xmx5g", "-jar", "/SAPP-2.0.jar", "-kofamscan", ]

arguments:
  - prefix: "-output"
    valueFrom: $(inputs.identifier).kofamscan.ttl
  - prefix: "-profile"
    valueFrom: /profiles/prokaryote.hal
  - prefix: "-kolist"
    valueFrom: /ko_list

outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).kofamscan.ttl

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