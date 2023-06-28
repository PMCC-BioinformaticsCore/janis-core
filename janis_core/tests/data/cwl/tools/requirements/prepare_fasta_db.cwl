#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Prepare fasta DB"

doc: |
    Prepares fasta file for so it does not contain duplicate fasta headers.
    Only looks at the first part of the header before any whitespace.
    Adds and incremental number in the header.

    Expects fasta file(s) or plaintext fasta(s). Not mixed!    

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "prepare_fasta_db"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          echo -e "\
          #/usr/bin/python3
          import sys\n\
          headers = set()\n\
          c = 0\n\
          for line in sys.stdin:\n\
            splitline = line.split()\n\
            if line[0] == '>':    \n\
              if splitline[0] in headers:\n\
                c += 1\n\
                print(splitline[0]+'.x'+str(c)+' '+' '.join(splitline[1:]))\n\
              else:\n\
                print(line.strip())\n\
              headers.add(splitline[0])\n\
            else:\n\
              print(line.strip())" > ./dup.py
          out_name=$1
          shift

          if file $@ | grep gzip; then
            zcat $@ | python3 ./dup.py | gzip > $out_name
          else
            cat $@ | python3 ./dup.py | gzip > $out_name
          fi

hints:
  SoftwareRequirement:
    packages:
      python3:
        version: ["3.10.6"]
        specs: ["https://anaconda.org/conda-forge/python"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/python:3.10.6
    
baseCommand: ["bash", "script.sh"]

inputs:
  fasta_files:
    type: File[]?
    label: fasta files
    doc: Fasta file(s) to be the prepared. Can also be gzipped (not mixe)
    inputBinding:
      position: 2
  output_file_name:
    type: string
    label: Output outfile
    inputBinding:
      position: 1

outputs:
  fasta_db:
    type: File?
    outputBinding:
      glob: $(inputs.output_file_name)

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
s:dateCreated: "2022-07-00"
s:dateModified: "2023-01-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/