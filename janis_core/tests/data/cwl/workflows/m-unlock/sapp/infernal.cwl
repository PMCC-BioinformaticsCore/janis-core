#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Infernal rna prediction"

doc: |
    Runs microbial infernal rna prediction on GBOL RDF file

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
  clanin:
    type: string
    doc: Path to Rfam.clanin file
    label: clanin file
    default: /unlock/references/databases/infernal/Rfam.clanin
    inputBinding:
      prefix: -clanin
  cm:
    type: string
    doc: Path to Rfam.cm file that is already indexed with cmpress
    label: cm file
    default: /unlock/references/databases/infernal/Rfam.cm
    inputBinding:
      prefix: -cm
  cmscan:
    type: string
    doc: Path to cmscan binary
    label: cmscan path
    default: /unlock/infrastructure/binaries/infernal/infernal-1.1.4-linux-intel-gcc/binaries/cmscan
    inputBinding:
      prefix: -cmscan
  threads:
    type: int?
    doc: number of threads to use for computational processes
    label: number of threads
    default: 1
    inputBinding:
      prefix: -cpu

baseCommand: ["java", "-Xmx5g", "-jar", "/SAPP-2.0.jar", "-infernal"]

arguments:
  - prefix: "-output"
    valueFrom: $(inputs.identifier).infernal.ttl

outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).infernal.ttl


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
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
  