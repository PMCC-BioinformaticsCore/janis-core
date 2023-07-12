#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Pairwise structural alignemnt using Kpax
doc: |
  Align crossmapped averaged structures against core average domain structure pairwise using Kpax
  Outputs a csv file with all the scores from pairwise alignments

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.result_file === 'string') {
          return [
            {"class": "File", "basename": inputs.result_file, "contents": "", writable: true}]; } 
        else { return [ inputs.result_file] ; } 
       }
    - entryname: pairwise_aligner.py
      entry:
        $include: Python/pairwise_aligner.py
      writable: false

hints:   
  SoftwareRequirement:
    packages:
      kpax:
        specs: [ "http://kpax.loria.fr/" ]      #link for the resgistration page of this tool is recommended (bio.tools?)
        version: [ "5.1.3.x64" ]


inputs:
  cath_fam_avg:
    label: CATH domain StIs to align 
    type:
      - "null"
      - File[]
      - File
    inputBinding:
      position: 1
      prefix: -c

  pfam_fam_avg:
    label: Pfam domain StIs to align 
    type:
      - "null"
      - File[]
      - File
    inputBinding:
      position: 2
      prefix: -p

  query_dir:
    type: Directory? 
    inputBinding:
      position: 3
      prefix: -d

  core_avg:
    label: Core average structure 
    type:
      - type: array
        items: File
      - File
    format: edam:format_1476
    inputBinding:
      position: 4
      itemSeparator: " "
      prefix: -t

  result_file:
    type: [ File?, string?]
    label: The output file
    default: "align_Struct_analysis.csv" 
    inputBinding:
      position: 5
      prefix: -r

  iteration:
    type: int
    default: 0

outputs:
  alignment_out:
    label: Alignment results from Kpax
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.result_file === 'string') {return inputs.result_file} else {return [ inputs.result_file.basename]}}


baseCommand:
  - python3
  - pairwise_aligner.py

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-7025-2241
    s:email: mailto:hbdhondge@gmail.com
    s:name: Hrishikesh Dhondge

  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-7035-3042
    s:email: mailto:isaure.chauvot-de-beauchene@loria.fr
    s:name: Isaure Chauvot de Beauchêne

  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-0399-8713
    s:email: mailto:marie-dominique.devignes@loria.fr
    s:name: Marie-Dominique Devignes

# s:citation: doi
s:codeRepository: https://gitlab.inria.fr/capsid.public_codes/CroMaSt
s:dateCreated: "2022-08-01"
s:license: https://mit-license.org/

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schemaorg-current-https.rdf
 - http://edamontology.org/EDAM_1.18.owl
