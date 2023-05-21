#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Download all the source files required by CroMaSt
doc: |
  This tool will download all the required files and store in the Data directory (default).


requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  MultipleInputFeatureRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: run_download.sh
      entry: |-
        mkdir $(inputs.data_directory)
        wget -O $(inputs.data_directory)/$(inputs.pfam_name).gz  $(inputs.pfam_url)
        gunzip $(inputs.data_directory)/$(inputs.pfam_name).gz
        wget -O $(inputs.data_directory)/$(inputs.cath_name)  $(inputs.cath_url)
        wget -O $(inputs.data_directory)/$(inputs.pdb_obs_name)  $(inputs.pdb_obs_url)


inputs:
  data_directory:
    label: Name for the Data directory
    type: string
    default: Data

  pfam_url:
    label: URL to download the file from
    type: string
  
  pfam_name:
    label: Name to save the file
    type: string
    default: pdbmap

  cath_url:
    label: URL to download CATH source file 
    type: string

  cath_name:
    label: Name for CATH source file
    type: string?
    default: cath-domain-description-file.txt

  pdb_obs_url:
    label: URL to download list of obsolete PDB entries 
    type: string

  pdb_obs_name:
    label: Name for list of obsolete PDB entries
    type: string?
    default: obsolete_PDB_entry_ids.txt


outputs:
  data_dir:
    type: Directory
    label: Data directory with all input data 
    outputBinding:
      glob: Data


baseCommand:
  - bash
  - run_download.sh 


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
