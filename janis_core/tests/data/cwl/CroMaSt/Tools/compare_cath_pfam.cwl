#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Compare residue-mapped instances from Pfam and CATH
doc: |
  Find the intersection between residue-mapped instances of Pfam and CATH lists.
  Allows variable domain boundaries in a certain range +/- 30aa. Produces three files: common domain instances,
  and unique domain instances to each Pfam and CATH.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.truedomains_file === 'string') {
          return [
            {"class": "File", "basename": inputs.truedomains_file, "contents": "{}", writable: true}]; } 
        else { return [ inputs.truedomains_file] ; } 
       }
    - entryname: compare_cath_pfam.py
      entry:
        $include: Python/compare_cath_pfam.py
      writable: false

inputs:
  resmapped_cath:
    label: All residue-mapped domain StIs with domain labels from CATH
    type: File
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -c

  resmapped_pfam:
    label: All residue-mapped domain StIs with domain labels from Pfam
    type: File
    format: edam:format_3752
    inputBinding:
      position: 2
      prefix: -p

  truedomains_file:
    type: [ File, string, "null"]
    label: Filename to store True domain StIs (.json)
    inputBinding:
      position: 3
      prefix: -f

  unique_pfam_struct:
    type: string?
    label: Filename for unique domain StIs from Pfam 
    default: unique_pfam.csv
    inputBinding:
      position: 4
      prefix: -uq_pf

  unique_cath_struct:
    type: string?
    label: Filename for unique domain StIs from CATH 
    default: unique_cath.csv
    inputBinding:
      position: 5
      prefix: -uq_ca

  min_dom_len:
    type: int
    label: Minimum domain length criteria to filter StIs
    default: 31
    inputBinding:
      prefix: -l


outputs:
  common_domains:
    label: True domain StIs per iteration
    type: File
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.truedomains_file === 'string') {return inputs.truedomains_file} else {return [ inputs.truedomains_file.basename]}}

  pfam_unique:
    label: Pfam domain StIs that are not in list of CATH domain StIs
    type: File
    format: edam:format_3752
    outputBinding:
      glob: $(inputs.unique_pfam_struct)

  cath_unique:
    label: CATH domain StIs that are not in list of Pfam domain StIs
    type: File
    format: edam:format_3752
    outputBinding:
      glob: $(inputs.unique_cath_struct)


baseCommand:
- python3
- compare_cath_pfam.py


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
