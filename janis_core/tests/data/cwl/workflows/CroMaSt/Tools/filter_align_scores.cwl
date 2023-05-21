#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Filter all structural instances for given Pfam families. 
doc: |
  The tool filter raw files from Pfam to retrieves all the available structural instances from the given Pfam families. 
  cwl-runner --cachedir=tmp_files/ --outdir=Results/ Workflow/separate_structures.cwl yml/separate_structures.yml 


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.passed_structs === 'string') {
          return [
            {"class": "File", "basename": inputs.passed_structs, "contents": "", writable: true}]; } 
        else { return [ inputs.passed_structs] ; } 
       }
    - |
      ${ 
        if (typeof inputs.failed_structs === 'string') {
          return [
            {"class": "File", "basename": inputs.failed_structs, "contents": "", writable: true}]; } 
        else { return [ inputs.failed_structs] ; } 
       }
    - entryname: filter_align_scores.py
      entry:
        $include: Python/filter_align_scores.py
      writable: false


inputs:
  aln_result:
    type: File
    label: Kpax alignment result
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -i

  all_struct_list:
    type: File
    format: edam:format_3464
    label: Json file as list of StIs 
    inputBinding:
      position: 2
      prefix: -x

  passed_structs:
    type: [ File, string, "null" ]
    label: Filename for list of domain StIs passing the threshold
    default: passed_structures_list.csv
    inputBinding:
      prefix: -p
  
  failed_structs:
    type: [ File, string, "null" ]
    label: Filename for list of domain StIs failed to pass the threshold
    default: failed_structures_list.csv
    inputBinding:
      prefix: -f
  
  aln_score:
    type: string?
    label: Kpax score type to use for filtering structures 
    default: Mscore
    inputBinding:
      prefix: -s

  threshold_val:
    type: float?
    label: The threshold to use for given aln_score
    default: 0.6
    inputBinding:
      prefix: -t


outputs:
  passed_structs_list:
    label: List of all domain StIs passing the threshold for given score
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.passed_structs === 'string') {return inputs.passed_structs} else {return [ inputs.passed_structs.basename]}}
  
  failed_structs_list:
    label: List of all domain StIs failed to pass the threshold for given score
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.failed_structs === 'string') {return inputs.failed_structs} else {return [ inputs.failed_structs.basename]}}


baseCommand:
- python3
- filter_align_scores.py


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
