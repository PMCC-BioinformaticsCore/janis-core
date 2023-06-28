#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      nanoqc :
        version: ["0.9.4"]
        specs: ["https://anaconda.org/bioconda/nanoqc"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/nanoqc:0.9.4

# baseCommand: ["bash", "script.sh"]
baseCommand: ["nanoQC"]

label: "NanoQC for Quality Check of nanopore FASTQ reads"
doc: |
    NanoQC for Quality Check of nanopore FASTQ reads

requirements:
   InlineJavascriptRequirement: {}
  #  InitialWorkDirRequirement:
    #  listing:
      #  - entryname: script.sh
        #  entry: |-
           #!/bin/bash
          #  source /root/miniconda/bin/activate
          #  conda init bash
          #  conda activate /unlock/infrastructure/conda/nanoqc/nanoqc_v0.9.4
          #  /unlock/infrastructure/conda/nanoqc/nanoqc_v0.9.4/bin/nanoQC $@

arguments: 
  - valueFrom: $(inputs.fastq)
  - valueFrom: $(inputs.outdir)
    prefix: --outdir

inputs:
  fastq:
    type: File
    doc: reads data in fastq.gz format
  outdir:
    type: string
    doc: name of output directory to be created
    inputBinding:
      prefix: --outdir

outputs:
  nanoqc_out: 
    type: Directory
    outputBinding:
      glob: $(inputs.outdir)

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
s:dateCreated: "2022-02-06"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/