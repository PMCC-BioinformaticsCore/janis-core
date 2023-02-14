#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Genome conversion"

doc: |
    Runs Genome conversion tool from SAPP

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
  embl:
    type: File?
    doc: Reference genome file used in EMBL format
    label: Reference genome
    inputBinding:
      prefix: -input
  fasta:
    type: File?
    doc: Reference genome file used in fasta format
    label: Reference genome
    inputBinding:
      prefix: -input
  identifier:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
    inputBinding:
      prefix: -id
  codon:
    type: int?
    doc: Codon table used for this organism
    label: Codon table
    inputBinding:
      prefix: -codon


baseCommand: ["java", "-Xmx5g", "-jar", "/SAPP-2.0.jar"]
arguments:
  - valueFrom: |
      ${
        if (inputs.embl) {
          return '-embl2rdf';
        }
        if (inputs.fasta) {
          return '-fasta2rdf';
        }
      }
  - valueFrom: |
      ${
        if (inputs.fasta) {
          return '-genome';
        } 
        return null;
      }
  - prefix: "-output"
    valueFrom: $(inputs.identifier).ttl

outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).ttl

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
  
