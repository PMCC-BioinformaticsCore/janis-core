#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "InterProScan annotation"

doc: |
    Runs Genome annotation with InterProScan using an RDF genome file according to the GBOL ontology

hints:
  SoftwareRequirement:
    packages:
      sapp:
        version: ["17.0.3"]
        specs: ["https://anaconda.org/conda-forge/openjdk"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/sapp:2.0

requirements:
  ResourceRequirement:
    ramMax: 10000
    coresMax: 10
  InlineJavascriptRequirement: {}

inputs:
  input:
    type: File
    doc: Reference genome file used in RDF format
    label: Reference genome
    inputBinding:
      prefix: -input
  cpu:
    type: int
    doc: Number of CPUs acccessible for InterProScan
    label: cpu limit
    inputBinding:
      prefix: -cpu
    default: 2
  applications:
    type: string
    doc: Which applications to use (e.g. TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,ProSiteProfiles,SMART,PRINTS,ProSitePatterns,Pfam,ProDom,MobiDBLite,PIRSF)
    label: applications
    inputBinding:
      prefix: -a
    default: Pfam
  identifier:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
  interpro:
    type: Directory
    doc: Path to the interproscan.sh folder
    label: InterProScan folder
    inputBinding:
      prefix: -path

baseCommand: ["java", "-Xmx5g", "-jar", "/SAPP-2.0.jar", "-interpro"]

arguments:
  - prefix: "-output"
    valueFrom: $(inputs.identifier).hdt

outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.identifier).hdt

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
  