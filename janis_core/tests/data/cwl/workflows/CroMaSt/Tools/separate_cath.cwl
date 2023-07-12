#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Filter all structural instances for given CATH superfamilies
doc: |
  The tool filter raw files from CATH to retrieve all the available structural instances from the given CATH superfamilies. 
  cwl-runner --cachedir=tmp_files/ --outdir=Results/ Workflow/separate_structures.cwl yml/separate_structures.yml 


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - '$(inputs.obs_pdb)'
    - '$(inputs.raw_cath)'
    - |
      ${ 
        if (typeof inputs.separate_cath === 'string') {
          return [
            {"class": "File", "basename": inputs.separate_cath, "contents": "", writable: true}]; } 
        else { return [ inputs.separate_cath] ; } 
       }
    - |
      ${ 
        if (typeof inputs.obsolete_cath === 'string') {
          return [
            {"class": "File", "basename": inputs.obsolete_cath, "contents": "", writable: true}]; } 
        else { return [ inputs.obsolete_cath] ; } 
       }
    - entryname: separate_cath.py
      entry:
        $include: Python/separate_cath.py
      writable: false
    - entryname: get_family_ids.py
      entry:
        $include: Python/get_family_ids.py
      writable: false
    - entryname: separate_pfam.py
      entry:
        $include: Python/separate_pfam.py
      writable: false


inputs:
  track_fams:
    type: File
    label: Family IDs per iteration
    format: edam:format_3464
    inputBinding:
      position: 1
      prefix: -f

  obs_pdb:
    type: File?
    label: All obsolete (deleted) PDB IDs 
    default: 
      class: File
      location: '../Data/obsolete_PDB_entry_ids.txt'
      basename: obsolete_PDB_entry_ids.txt
    inputBinding:
      prefix: -d

  raw_cath:
    type: File?
    label: Raw file from CATH with all domain instances
    default: 
      class: File
      location: '../Data/cath-domain-description-file.txt'
      basename: cath-domain-description-file.txt
    inputBinding:
      prefix: -c

  separate_cath:
    type: [ File?, string?]
    label: Filename for filtered structures from CATH
    default: Filtered_CATH.csv
    inputBinding:
      prefix: -n
  
  obsolete_cath:
    type: [ File?, string?]
    label: Filename for obsolete cath structures 
    default: obsolete_cath.txt
    inputBinding:
      prefix: -o

  split_suffix:
    type: string?
    label: Suffix for splitted files (.csv) used to parallelize task
    default: part.csv
    inputBinding:
      prefix: -s

  min_dom_len:
    type: int
    label: Minimum domain length criteria to filter structural instances
    default: 31
    inputBinding:
      prefix: -l


outputs:
  cath_structs:
    type: File
    label:  A file containing all the filtered structures from CATH
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.separate_cath === 'string') {return inputs.separate_cath} else {return [ inputs.separate_cath.basename]}}

  cath_obs:
    type: File
    format: edam:format_2330
    label:  Obsolete structures from the list of given CATH families
    outputBinding:
      glob: ${ if (typeof inputs.obsolete_cath === 'string') {return inputs.obsolete_cath} else {return [ inputs.obsolete_cath.basename]}}

  splitted_cath_sep:
    type: File[]
    label: list of files containing all the filtered structures from CATH 
    format: edam:format_3752
    outputBinding:
      glob: '*$(inputs.split_suffix)'

baseCommand:
- python3
- separate_cath.py


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
