#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Checks the alignment score for given threshold
doc: |
  Checks the alignment score for each aligned structure based on the given threshold
  Outputs the structural instances passing and failing the threshold in separate files

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: check_threshold.py
      entry:
        $include: Python/check_threshold.py
      writable: false

inputs: 
  aln_res_file:
    label: Alignment result from Kpax
    type: File
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -a

  fam_tracker:
    label: Family Ids per iteration (.json)
    type: File
    format: edam:format_3464
    inputBinding:
      position: 2
      prefix: -f

  aln_score:
    label: Kpax score type to use for filtering structures
    type: string?
    default: Mscore
    inputBinding:
      position: 3
      prefix: -s

  threshold_val:
    label: The threshold to use for given aln_score
    type: float?
    default: 0.6
    inputBinding:
      position: 4
      prefix: -t

  pfam_crossmap:
    label: All Pfam domain StIs corresponding to cross-mapped family
    type: File?
    format: edam:format_3464
    inputBinding:
      position: 5
      prefix: -px
  
  cath_crossmap:
    label: All CATH domain StIs corresponding to cross-mapped family
    type: File?
    format: edam:format_3464
    inputBinding:
      position: 6
      prefix: -cx
  

outputs:
  pfam_crossmap_passed:
    label: Cross-mapped families with Pfam domain StIs passing the threshold
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: "crossmapped_pfam_passed.json"
  
  pfam_crossmap_failed:
    label: Cross-mapped families with Pfam domain StIs failed to pass the threshold
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: "crossmapped_pfam_failed.json"

  cath_crossmap_passed:
    label: Cross-mapped families with CATH domain StIs passing the threshold
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: "crossmapped_cath_passed.json"
  
  cath_crossmap_failed:
    label: Cross-mapped families with CATH domain StIs failed to pass the threshold
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: "crossmapped_cath_failed.json"


baseCommand:
  - python3
  - check_threshold.py

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
