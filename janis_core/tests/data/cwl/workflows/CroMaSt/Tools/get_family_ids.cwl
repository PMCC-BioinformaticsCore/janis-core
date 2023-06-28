#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Get domain family ids
doc: |
  Get domain family ids from CATH and Pfam databases from parameter file provided by user
  and write these ids to a separate file.


requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  MultipleInputFeatureRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.fam_tracker === 'string') {
          return [
            {"class": "File", "basename": inputs.fam_tracker, "contents": "{}", writable: true}]; } 
        else { return [ inputs.fam_tracker] ; } 
       }
    - entryname: get_family_ids.py
      entry:
        $include: Python/get_family_ids.py
      writable: false


inputs:
  pfam_ids:
    type: string[]?
    label: List of Pfam family IDs
    inputBinding:
      position: 1
      prefix: -p

  cath_ids:
    type: string[]?
    label: List of CATH family (node) IDs
    inputBinding:
      position: 2
      prefix: -c

  iteration_no:
    type: int
    label: Iteration number to keep track of the workflow
    inputBinding:
      position: 3
      prefix: -n

  fam_tracker:
    type: [ File, string, "null" ]
    label: File to track family IDs per iteration
    default: family_ids.json
    inputBinding:
      position: 4
      prefix: -f


outputs:
  family_ids:
    type: File
    label: File with information for family IDs per iteration
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.fam_tracker === 'string') {return inputs.fam_tracker} else {return [ inputs.fam_tracker.basename]}}


baseCommand:
  - python3
  - get_family_ids.py


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
