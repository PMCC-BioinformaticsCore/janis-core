#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: Produce a list of residue-mapped structural domain instances from CATH ids
doc: |
  Retrieve and process the PDB structures corresponding to the CATH superfamily ids resulting in a list of residue-mapped
  structural domain instances along with lost structural instances
  (requires Data/cath_domain_description_file.txt downloaded from CATH and uses SIFTS resource for PDB to UniProt residue Mapping)

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  family_idsfile:
    label: File with the family IDs per iteration 
    type: File
    format: edam:format_3464
  siftsdir:
    label: Directory for storing all SIFTS files
    type: Directory
  resmapped_file: 
    label: Filename for CATH inconsistent domain StIs
    type: string
  lost_merged: 
    label: Filename for Pfam inconsistent domain StIs
    type: string
  min_dom_size: 
    label: Threshold for minimum domain length
    type: int


steps:
  filter_cath_structures:
    label: Filter all domain StIs for given CATH superfamilies
    run: 'separate_cath.cwl'
    in:
      track_fams: family_idsfile
      min_dom_len: min_dom_size
    out: [cath_structs, cath_obs, splitted_cath_sep]

  resmapping_cath_structs:
    run:
      class: Workflow
      label: Mapping of residue numbering from PDB to UniProt
      inputs:
        flt_files: 
          type: File
          format: edam:format_3752
        sifts: 
          type: Directory

      steps:
        resmapping_for_CATH_PDB2UP:
          run: 'resmapping_cath2up.cwl'
          in:
            cath_sep: flt_files
            sifts_dir: sifts
          out: [cath_resmapped, cath_lost]

      outputs:
        resmapped_cath:
          type: File
          format: edam:format_3752
          outputSource: resmapping_for_CATH_PDB2UP/cath_resmapped
              
        lost_insta_cath:
          type: File
          format: edam:format_2330
          outputSource: resmapping_for_CATH_PDB2UP/cath_lost

    scatter: flt_files
    in:
      flt_files: filter_cath_structures/splitted_cath_sep
      sifts: siftsdir
    out: [ resmapped_cath, lost_insta_cath ]
        
  add_domain_positions:
    label: Add domain positions to residue-mapped StIs
    run: 'add_domain_num.cwl'
    in: 
      resmapped_files: resmapping_cath_structs/resmapped_cath
      dom_posi_file: resmapped_file
    when: $(inputs.resmapped_files.length > 0)
    out: [ resmapped_domains ]

  collect_lost_instances:
    label: Merge obsolete and inconsistent domain StIs together
    run: 'gather_lost_resmap.cwl'
    in:
      lost_instance: resmapping_cath_structs/lost_insta_cath
      obs_insta: filter_cath_structures/cath_obs
      outfile: lost_merged
    out: [ lost_domain_list ]


outputs: 
  cath_domain_posi_file: 
    label: All residue-mapped domain StIs with domain labels
    type: File
    format: edam:format_3752
    outputSource: add_domain_positions/resmapped_domains

  cath_total_lost_structures:
    label: Obsolete and inconsistent domain StIs together 
    type: File
    format: edam:format_3464
    outputSource: collect_lost_instances/lost_domain_list

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
