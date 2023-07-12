#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Add cross-mapped structural instances to residue-mapped structural instances
doc: |
  Add crossmapped domain instances from last iteration to current list of residue mapped domain instances.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.pfam_result === 'string') {
          return [
            {"class": "File", "basename": inputs.pfam_result, "contents": "", writable: true}]; } 
        else { return [ inputs.pfam_result] ; } 
       }
    - |
      ${ 
        if (typeof inputs.cath_result === 'string') {
          return [
            {"class": "File", "basename": inputs.cath_result, "contents": "", writable: true}]; } 
        else { return [ inputs.cath_result] ; } 
       }
    - entryname: add_crossmapped2resmapped.py
      entry:
        $include: Python/add_crossmapped2resmapped.py
      writable: false

inputs:
  pfam_resmapped:
    type: File?
    label: Pfam residue-mapped domain StIs from current iteration
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -p

  cath_resmapped:
    type: File?
    label: CATH residue-mapped domain StIs from current iteration
    format: edam:format_3752
    inputBinding:
      position: 2
      prefix: -c

  pfam_crossmapped:
    type: File?
    label: Pfam cross-mapped domain StIs from previous iteration
    format: edam:format_3464
    inputBinding:
      position: 3
      prefix: -px

  cath_crossmapped:
    type: File?
    label: CATH cross-mapped domain StIs from previous iteration
    format: edam:format_3464
    inputBinding:
      position: 4
      prefix: -cx

  pfam_result:
    type: [ File?, string?]
    label: Filename to merge cross-mapped and residue-mapped domain StIs from Pfam
    default: 'pfam_res_crossMapped.csv'
    inputBinding:
      position: 5
      prefix: -pr

  cath_result:
    type: [ File?, string?]
    label: Filename to merge cross-mapped and residue-mapped domain StIs from CATH
    default: 'cath_res_crossMapped.csv'
    inputBinding:
      position: 6
      prefix: -cr


outputs:
  pfam_structs:
    label: Merged cross-mapped and residue-mapped domain StIs from Pfam
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.pfam_result === 'string') {return inputs.pfam_result} else {return inputs.pfam_result.basename}}

  cath_structs:
    label: Merged cross-mapped and residue-mapped domain StIs from CATH
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.cath_result === 'string') {return inputs.cath_result} else {return inputs.cath_result.basename}}


baseCommand:
  - python3
  - add_crossmapped2resmapped.py

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
