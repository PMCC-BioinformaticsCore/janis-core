#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      krona :
        version: ["2.8.1"]
        specs: ["https://anaconda.org/bioconda/krona"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/krona:2.8.1

baseCommand: [ ktImportTaxonomy ]

label: Krona
doc: |
    Visualization of Kraken2 report results.
    ktImportText -o $1 $2

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "krona_output"
        writable: true
      # - entryname: script.sh
        # entry: |-
          #!/bin/bash
          # source /root/miniconda/bin/activate
          # conda init bash
          # conda activate /unlock/infrastructure/conda/krona/krona_v2.8.1
          # ktImportTaxonomy -t 5 -m 3 $@

arguments:
  - prefix: "-o"
    valueFrom: krona_output/$(inputs.kraken.nameroot)_krona.html

inputs:
  kraken:
    type: File
    label: Tab-delimited text file
    inputBinding:
      position: 10 # Krona requires the input after the ouput flag
  taxonomy:
    type: int
    label: Taxon column
    doc: Column number for taxon information (default for kraken)
    default: 5
    inputBinding:
      position: 1
      prefix: -t
  counts:
    type: int
    label: Counts column
    doc: Column number for count information (default for kraken)
    default: 3
    inputBinding:
      position: 2
      prefix: -m

outputs:
  # krona_outdir:
  #   type: Directory
  #   outputBinding:
  #     glob: krona_output
  krona_html:
    type: File
    outputBinding:
      glob: krona_output/$(inputs.kraken.nameroot)_krona.html

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
s:dateCreated: "2021-12-10"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
