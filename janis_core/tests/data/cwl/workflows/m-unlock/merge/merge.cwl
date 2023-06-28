#!/usr/bin/env cwltool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ cat ]

label: "FASTQ merge tool"
doc: |
    Performs a cat on multiple sets of files

requirements:
 - class: ShellCommandRequirement
 - class: InlineJavascriptRequirement

arguments:
  - shellQuote: false
    valueFrom: >
      ${
        var cmd = inputs.files.path + "/*.fastq > " + inputs.identifier + ".fastq";
        return cmd;
      }

inputs:
  identifier:
    type: string
    doc: Name of the output file
  files:
    type: Directory
    doc: folder with files to concatenate

outputs:
  merged_fastq:
    type: stdout
  # fastq_outdir: #UNCOMMENT FOR SPECIFIC OUTPUT FOLDER
  #   type: Directory
  #   outputBinding:
  #     glob: .
  #     outputEval: |
  #       ${
  #         self[0].basename = "fastq_merged";
  #         return self[0]
  #       }

stdout: $(inputs.identifier).fastq

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
s:dateCreated: "2020-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
  