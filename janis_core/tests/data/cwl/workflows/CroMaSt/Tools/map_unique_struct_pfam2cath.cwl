#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Map unique Pfam instances to CATH db
doc: |
  Maps the unique instances from Pfam to the whole CATH database 
  (using residue numbering from PDB allowing variable domain boundaries, by default +/-30aa)

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.crossmap_pfam === 'string') {
          return [
            {"class": "File", "basename": inputs.crossmap_pfam, "contents": "", writable: true}]; } 
        else { return [ inputs.crossmap_pfam] ; } 
       }
    - |
      ${ 
        if (typeof inputs.no_crossmap === 'string') {
          return [
            {"class": "File", "basename": inputs.no_crossmap, "contents": "", writable: true}]; } 
        else { return [ inputs.no_crossmap] ; } 
       }
    - entryname: map_unique_struct_pfam2cath.py
      entry:
        $include: Python/map_unique_struct_pfam2cath.py
      writable: false

inputs:
  pfam_unq:
    label: Pfam domain StIs to cross-map against whole CATH db
    type: File
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -p

  cath_raw:
    type: File?
    label: Raw file from CATH with all domain instances
    default: 
      class: File
      location: '../Data/cath-domain-description-file.txt'
      basename: cath-domain-description-file.txt
    inputBinding:
      prefix: -c

  crossmap_pfam:
    type: [ File?, string?]
    label: User-defined filename for filtered structures from CATH
    default: pfam_crossMapped_cath.jsonx
    inputBinding:
      prefix: -x
  
  no_crossmap:
    type: [ File?, string?]
    label: User-defined filename for unampped domain StIs 
    default: pfam_unq_unmapped.jsonx
    inputBinding:
      prefix: -u

  min_dom_len:
    type: int
    label: Minimum domain length criteria to filter domain StIs
    default: 31
    inputBinding:
      prefix: -l


outputs:
  pfam_crossmapped:
    label: Pfam cross-mapped domin StIs family-wise
    type:
      type: array
      items: File
    format: edam:format_3464
    outputBinding:
      glob: "*.json"

  allcrossmap_pfam:
    label: All Pfam cross-mapped domin StIs family-wise together
    type: File
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.crossmap_pfam === 'string') {return inputs.crossmap_pfam} else {return [ inputs.crossmap_pfam.basename]}}

  pfam_unmapped:
    label: All Pfam un-mapped domin StIs 
    type: File
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.no_crossmap === 'string') {return inputs.no_crossmap} else {return [ inputs.no_crossmap.basename]}}

baseCommand:
- python3
- map_unique_struct_pfam2cath.py

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
