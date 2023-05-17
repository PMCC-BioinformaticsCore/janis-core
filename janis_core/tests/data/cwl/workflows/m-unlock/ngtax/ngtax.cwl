#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "NGTax amplicon analysis"

doc: |
    Runs NGTAX amplicon analysis

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference_db)
hints:
  SoftwareRequirement:
    packages:
      java:
        version: ["17.0.3"]
        specs: ["https://anaconda.org/conda-forge/openjdk"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/ngtax:2.2.9

inputs:
  forward_primer: # "[AG]GGATTAGATACCC"
    type: string
    doc: Forward primer used
    label: The forward primer used
    inputBinding:
      prefix: -for_p
  reverse_primer: # "CGAC[AG][AG]CCATGCA[ACGT]CACCT"
    type: string?
    doc: Reverse primer used
    label: The reverse primer used
    inputBinding:
      prefix: -rev_p
  reference_db:
    type: File?
    doc: Reference database used in FASTA format
    label: Reference database
    inputBinding:
      prefix: -refdb
  folder:
    type: Directory?
    doc: Folder containing demultiplexed files_to_folder
    label: Demultiplexed folder
    inputBinding:
      prefix: -folder
  rev_read_len: 
    type: int?
    doc: Read length of the reverse read
    label: Reverse read length
    inputBinding:
      prefix: -rev_read_len
  for_read_len: 
    type: int
    doc: Read length of the reverse read
    label: Reverse read length
    inputBinding:
      prefix: -for_read_len
  # minimum_threshold: 
    # type: float
    # doc: Minimum threshold detectable, expressed in percentage
    # label: Minimum threshold
    # inputBinding:
      # prefix: -minimumThreshold
  mock3:
    type: string?
    doc: Mock3 reference selection
    label: Mock3 reference
  mock4:
    type: string?
    doc: Mock4 reference selection
    label: Mock4 reference
  sample:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
  fragment:
    type: string
    doc: Subfragment that is being analysed (e.g. V1-V3 or V5-region)
    label: Subfragment name
  primersRemoved:
    type: boolean?
    doc: Wether the primers are removed or not from the input files
    label: Primers are removed
    inputBinding:
      prefix: -primersRemoved

baseCommand: ["java", "-jar", "/NGTax-2.2.9.jar", "-ngtax", "-mapFile", "cwl_mapping_file.txt"]

arguments:
  - prefix: "-t"
    valueFrom: $(inputs.sample)_NG-Tax_$(inputs.for_read_len).ttl
  - prefix: "-b"
    valueFrom: $(inputs.sample)_NG-Tax_$(inputs.for_read_len).biom
  - prefix: "-mock3"
    valueFrom: $(inputs.mock3)
  - prefix: "-mock4"
    valueFrom: $(inputs.mock4)
  - prefix: "-fragment"
    valueFrom: $(inputs.fragment)

stdout: ngtax2.stdout.log

outputs:
  biom: 
    type: File
    outputBinding:
      glob: "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).biom"
  turtle: 
    type: File
    outputBinding:
      glob: "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).ttl"
  stdout_out:
    type: File
    outputBinding:
      glob: ngtax2.stdout.log

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:dateModified: "2023-02-06"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
