#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

# baseCommand: ["bash", "script.sh"]
baseCommand: [ metaquast.py ]

label: "metaQUAST: Quality Assessment Tool for Metagenome Assemblies"
doc: |
  Runs the Quality Assessment Tool for Metagenome Assemblies application
  
  Necessary to install the pre-release to prevent issues:
  https://github.com/ablab/quast/releases/tag/quast_5.1.0rc1
  
  The working installation followed the method in http://quast.sourceforge.net/docs/manual.html:
  $ wget https://github.com/ablab/quast/releases/download/quast_5.1.0rc1/quast-5.1.0rc1.tar.gz
  $ tar -xzf quast-5.1.0rc1.tar.gz
  $ cd quast-5.1.0rc1/
  $ ./setup.py install_full

hints:
  SoftwareRequirement:
    packages:
      quast:
        version: ["5.2.0"]
        specs: ["https://anaconda.org/bioconda/quast"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/quast:5.2.0

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "metaQUAST_results"
      writable: true
    # - entryname: script.sh
    #   entry: |-
    #     #!/bin/bash
    #     # source /root/miniconda/bin/activate
    #     # conda init bash
    #     # conda activate /unlock/infrastructure/conda/quast/quast_v5.2.0
    #     metaquast.py $@

arguments:
  - valueFrom: $(inputs.assembly)
  - valueFrom: "metaQUAST_results"
    prefix: --output-dir

    # If no custom BLAST database is given it will use the SILVA database pre-build in the docker container.
  - |
    ${
      if (inputs.blastdb){
        return ["--blast-db", inputs.blastdb.path];
      } else {
        return ["--blast-db", "/venv/silva_138.1/"];
      }
    }

inputs:
  assembly:
    type: File
    doc: The input assembly in fasta format
#   metaquast_out:
#     type: string
#     doc: pathway to output directory including name of folder
#     inputBinding:
#       prefix: --output-dir
  threads:
    type: int
    default: 1
    inputBinding:
      prefix: --threads
  blastdb:
    type: Directory?
    label: Reference BLASTdb
    doc: Reference BLASTdb. The path should point to a directory containing .nsq file or to .nsq file itself. (default is silva)
    inputBinding:
      prefix: --blast-db 
  space_efficient:
    type: boolean
    label: Run space efficient
    doc: Uses less diskspace. Removes aux files with some details for advanced analysis. (default true)
    inputBinding:
      prefix: --space-efficient
    default: true

outputs:
  metaquast_outdir:
    type: Directory
    outputBinding:
      glob: metaQUAST_results
      #glob: $(inputs.metaquast_out)
  meta_combined_ref:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/combined_reference
      #glob: $(inputs.metaquast_out)/combined_reference
  meta_icarusDir:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/icarus_viewers
  metaquast_krona:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/krona_charts
  not_aligned:  
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/not_aligned
  meta_downloaded_ref:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/quast_downloaded_references
  runs_per_reference:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/runs_per_reference
  meta_summary:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/summary
  meta_icarus:
    type: File?
    outputBinding:
      glob: metaQUAST_results/icarus.html
  metaquast_log:
    type: File?
    outputBinding:
      glob: metaQUAST_results/metaquast.log
  metaquast_report:
    type: File?
    outputBinding:
      glob: metaQUAST_results/report.html

# QUAST outputs, in case the tool does a normal QUAST instead of metaQUAST.
# if metaQUAST doesn't find a reference, it changes to QUAST.
  basicStats:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/basic_stats
  quast_icarusDir:
    type: Directory?
    outputBinding:
      glob: metaQUAST_results/icarus_viewers
  quast_icarusHtml:
    type: File?
    outputBinding:
      glob: metaQUAST_results/icarus.html
  quastReport:
    type: File[]?
    outputBinding:
      glob: metaQUAST_results/report.*
  quastLog:
    type: File?
    outputBinding:
      glob: metaQUAST_results/quast.log
  transposedReport:
    type: File[]?
    outputBinding:
      glob: metaQUAST_results/transposed_report.*

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
s:dateCreated: "2021-02-09"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
