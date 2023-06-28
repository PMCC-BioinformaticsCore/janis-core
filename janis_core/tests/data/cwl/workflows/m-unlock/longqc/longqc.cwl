#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      longqc :
        version: ["1.2.0"]
        specs: ["https://anaconda.org/bioconda/longqc"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/longqc:1.2.0

baseCommand: ["bash", "script.sh"]

label: "LongQC for Quality Check of nanopore FASTQ reads"
doc: |
    LongQC for Quality Check of nanopore FASTQ reads

    ***README: Unsolved issue with parsing file:
    https://github.com/yfukasawa/LongQC/issues/37
    "pandas.errors.EmptyDataError: No columns to parse from file"

requirements:
   InlineJavascriptRequirement: {}
   InitialWorkDirRequirement:
     listing:
       - entryname: script.sh
         entry: |-
           #!/bin/bash
           source /root/miniconda/bin/activate
           conda init bash
           conda activate /unlock/infrastructure/conda/longqc
           python3 /unlock/infrastructure/conda/longqc/bin/longQC.py sampleqc $@

arguments: 
  - valueFrom: ont-rapid
    prefix: --preset
  - valueFrom: $(inputs.sampleqc)

inputs:
  threads:
    type: int
    default: 4 #or it complains about needing 4 or higher
    inputBinding:
      prefix: --ncpu
  output_dir:
    type: Directory
    doc: output directory to be created with full path
    inputBinding:
      prefix: --output
  sampleqc:
    type: File
    doc: FASTQ input file

outputs:
  longqc_out: 
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-5516-8391
    s:email: mailto:german.royvalgarcia@wur.nl
    s:name: Germán Royval
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
s:dateCreated: "2022-01-31"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
