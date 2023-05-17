#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "SAPP conversion RDF2GTF"

doc: |
    SAPP conversion tool utilizing the function RDF2GTF

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
      prefix: "-input"

arguments:
  - prefix: "-output"
    valueFrom: $(inputs.inputFile.nameroot)

baseCommand: ["java", "-Xmx5G", "-jar", "/SAPP-2.0.jar", "-rdf2gtf"]

outputs:
  genomefasta:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.nameroot).fasta
  gtf:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.nameroot).gtf

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
s:dateModified: "2022-05-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
  