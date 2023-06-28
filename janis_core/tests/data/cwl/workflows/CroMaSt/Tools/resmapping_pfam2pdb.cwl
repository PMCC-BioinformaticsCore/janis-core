#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Residue mapping from UniPort to PDB numbering for Pfam structural instances
doc: |
  The tool maps the residues from UniPort to PDB numbering for all the Pfam structural instances using SIFTS mapping.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.resmapping_file === 'string') {
          return [
            {"class": "File", "basename": inputs.resmapping_file, "contents": "", writable: true}]; } 
        else { return [ inputs.resmapping_file] ; } 
       }
    - |
      ${ 
        if (typeof inputs.reslost === 'string') {
          return [
            {"class": "File", "basename": inputs.reslost, "contents": "", writable: true}]; } 
        else { return [ inputs.reslost] ; } 
       }
    - '$({class: "Directory", basename: inputs.sifts_dir.location, listing: []})'
    - entryname: resmapping_pfam2pdb.py
      entry:
        $include: Python/resmapping_pfam2pdb.py
      writable: false
    - entryname: residue_mapping.py
      entry:
        $include: Python/residue_mapping.py

inputs:
  pfam_sep:
    label: Filtered structures from Pfam
    type: File
    format: edam:format_3752
    inputBinding:
      position: 1
      prefix: -f

  sifts_dir:
    type: Directory
    inputBinding:
      position: 2
      prefix: -s
  
  resmapping_file:
    type: [ File?, string?]
    label: Filename for the residue mapped file
    default: pfam_resMapped.csv
    inputBinding:
      position: 3
      prefix: -m

  reslost:
    type: [ File?, string?]
    label: Filename for the lost structures while residue mapping
    default: lost_pfam.txt
    inputBinding:
      position: 4
      prefix: -l


outputs:
  pfam_resmapped:
    type: File
    format: edam:format_3752
    outputBinding:
      glob: ${ if (typeof inputs.resmapping_file === 'string') {return inputs.resmapping_file} else {return [ inputs.resmapping_file.basename]}}

  pfam_lost:
    type: File
    format: edam:format_2330
    outputBinding:
      glob: ${ if (typeof inputs.reslost === 'string') {return inputs.reslost} else {return [ inputs.reslost.basename]}}

baseCommand:
- python3
- resmapping_pfam2pdb.py


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
