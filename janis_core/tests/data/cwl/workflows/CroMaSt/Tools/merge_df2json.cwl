#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Merge all passed domain instances together 
doc: |
  Helps to collect all domain instances passed to this tool and outputs the list with all these domain instances together.

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  MultipleInputFeatureRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.unmapped_out === 'string') {
          return [
            {"class": "File", "basename": inputs.unmapped_out, "contents": "{}", writable: true}]; } 
        else { return [ inputs.unmapped_out] ; } 
       }
    - entryname: merge_unmapped.py
      entry:
        $include: Python/merge_df2json.py
      writable: false

inputs:
  pfam_unmapped:
    label: Un-mapped domain StIs from Pfam
    type: File?
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -p

  cath_unmapped:
    label: Un-mapped domain StIs from CATH
    type: File?
    format: edam:format_3752
    inputBinding:
      position: 2
      prefix: -c

  pfam_crossmapped:
    label: Cross-mapped domain StIs from Pfam only for failed StIs
    type: File?
    format: edam:format_3464
    inputBinding:
      position: 3
      prefix: -px

  cath_crossmapped:
    label: Cross-mapped domain StIs from CATH only for failed StIs
    type: File?
    format: edam:format_3464
    inputBinding:
      position: 4
      prefix: -cx

  unmapped_out:
    label: Output filename for either Domain-like or failed domain StIs
    type: [ File, string, "null" ]
    default: domain_StIs_f.json
    inputBinding:
      position: 3
      prefix: -o


outputs:
  unmapped_list:
    label: With total unmapped (passed or failed) domain StIs until this iteration
    type: File
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.unmapped_out === 'string') {return inputs.unmapped_out} else {return [ inputs.unmapped_out.basename]}}


baseCommand:
  - python3
  - merge_unmapped.py


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
