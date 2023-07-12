#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

hints:
  SoftwareRequirement:
    packages:
      kraken2:
        version: ["2.1.2"]
        specs: ["https://anaconda.org/bioconda/kraken2"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/kraken2:2.1.2

# baseCommand: [ /unlock/infrastructure/binaries/kraken2-2.0.9-beta/kraken2 ]
# baseCommand: ["bash","-x", "script.sh"]
baseCommand: [kraken2]

label: "Kraken2 metagenomics read classification"
doc: |
    Kraken2 metagenomics read classification.
    
    Updated databases available at: https://benlangmead.github.io/aws-indexes/k2 (e.g. PlusPF-8)
    Original db: https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: "$({class: 'Directory', listing: []})"
        entryname: "kraken2_output"
        writable: true
      - entryname: script.sh
        entry: |-
          #!/bin/bash
          PATH=/unlock/infrastructure/binaries/BLAST/ncbi-blast-2.13.0+/bin/:"$PATH"
          # chmod +x /unlock/infrastructure/binaries/kraken2/kraken2-v2.1.2/*
          DBNAME=KRAKEN2_$1
          threads=$2
          # /unlock/infrastructure/binaries/kraken2/kraken2-v2.1.2/
          kraken2-build --download-taxonomy --db $DBNAME
          shift;shift;
          for file in "$@"
          do
              echo "Processing $file";
              # /unlock/infrastructure/binaries/kraken2/kraken2-v2.1.2/
              kraken2-build --add-to-library $file --db $DBNAME
          done
          # /unlock/infrastructure/binaries/kraken2/kraken2-v2.1.2/
          kraken2-build --threads $threads --build --db $DBNAME


inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
    inputBinding:
      position: 1
  threads:
    type: int?
    default: 1
    inputBinding:
      position: 2
  references:
    type: string[]
    doc: List of genome references for kraken
    inputBinding:
      position: 10

outputs:
  kraken2_database: 
    type: Directory
    outputBinding:
      glob: KRAKEN2_*


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
s:dateCreated: "2022-01-01"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/