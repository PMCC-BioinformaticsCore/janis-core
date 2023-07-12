#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Create parameter file for next iteration
doc: |
  Create parameter file for next iteration from previous parameter file
  Filter the pairwise alignments to retrieve family ids passing the threshold for a given Kpax score type

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: runtime_parm.json
      entry: |-
        {
        "filename": $(inputs.fam_tracker),
        "true_domain_file": $(inputs.true_domains),
        "core_domain_struct": $(inputs.core_domain_struct),
        "crossmap_pfam": $(inputs.crossmap_pfam),
        "crossmap_cath": $(inputs.crossmap_cath),
        "domain_like": $(inputs.domain_like),
        "failed_domain": $(inputs.failed_domains),
        "parm_file": $(inputs.in_paramfile),
        "pfam_lost_structs": $(inputs.pfam_lost),
        "cath_lost_structs": $(inputs.cath_lost)
        }
    - |
      ${ 
        if (typeof inputs.next_paramfile === 'string') {
          return [
            {"class": "File", "basename": inputs.next_paramfile, "contents": "", writable: true}]; } 
        else { return [ inputs.next_paramfile] ; } 
       }
    - entryname: create_param.py
      entry:
        $include: Python/create_param.py
      writable: false

inputs:
  in_paramfile:
    label: Parameter file for current iteration
    type: File
    format: edam:format_3750
    inputBinding:
      position: 1
      prefix: -i

  next_paramfile:
    label: Parameter file for next iteration
    type: [ File?, string?]
    default: new_param.yml
    inputBinding:
      position: 2
      prefix: -o
  
  crossmap_pfam:
    label: cross-mapped families for Pfam domain StIs that passes the threshold
    type: File
    format: edam:format_3464
    inputBinding:
      position: 3
      prefix: -px
  
  crossmap_cath:
    label: cross-mapped families for CATH domain StIs that passes the threshold
    type: File
    format: edam:format_3464 
    inputBinding:
      position: 4
      prefix: -cx

  fam_tracker:
    label: Family ids per iteration
    type: File
    format: edam:format_3464
  
  true_domains:
    label: True domain StIs
    type: File
    format: edam:format_3464

  core_domain_struct:
    label: Core average structure
    type: 
      - File
      - File[] 
    format: edam:format_1476
  
  domain_like:
    label: Domain-like StIs
    type: File
    format: edam:format_3464

  failed_domains:
    label: Failed domain StIs
    type: File
    format: edam:format_3464

  pfam_lost:
    label: Obsolete and inconsistent Pfam domain StIs 
    type: [ File, "null"]

  cath_lost:
    label: Obsolete and inconsistent CATH domain StIs 
    type: [ File, "null"]

outputs:
  next_parmfile:
    label: Parameter file for next iteration
    type: File
    format: edam:format_3750
    outputBinding:
      glob: ${ if (typeof inputs.next_paramfile === 'string') {return inputs.next_paramfile} else {return inputs.next_paramfile.basename}}

baseCommand:
  - python3
  - create_param.py

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
