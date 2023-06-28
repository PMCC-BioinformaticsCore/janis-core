#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "SAPP conversion RDF2FASTA"

doc: |
    SAPP conversion tool utilizing the function RDF2FASTA

hints:
  SoftwareRequirement:
    packages:
      sapp:
        version: ["17.0.3"]
        specs: ["https://anaconda.org/conda-forge/openjdk"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/sapp:2.0

requirements:
 - class: InlineJavascriptRequirement

inputs:
  inputFile:
    type: File
    inputBinding:
      prefix: "-i"

arguments:
  - prefix: "-transcript"
    valueFrom: $(inputs.inputFile.nameroot)_transcripts.fasta
  - prefix: "-protein"
    valueFrom: $(inputs.inputFile.nameroot)_proteins.fasta

baseCommand: ["java", "-Xmx5G", "-jar", "/SAPP-2.0.jar", "-rdf2fasta"]

outputs:
  transcripts:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.nameroot)_transcripts.fasta
  proteins:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.nameroot)_proteins.fasta

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
s:dateCreated: "2020-08-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/