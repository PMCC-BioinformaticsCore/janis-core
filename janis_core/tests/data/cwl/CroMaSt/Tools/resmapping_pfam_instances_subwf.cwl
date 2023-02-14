#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: Produce a list of residue-mapped structural domain instances from Pfam ids
doc: |
  Retrieve and process the PDB structures corresponding to the Pfam family ids resulting in a list of residue-mapped
  structural domain instances along with lost structural instances
  (requires Data/pdbmap downloaded from Pfam and uses SIFTS resource for UniProt to PDB residue Mapping)

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
    label: Filename for CATH inconsistent structural instances
    type: string
  lost_merged: 
    label: Filename for Pfam inconsistent structural instances
    type: string
  min_dom_size: 
    label: Threshold for minimum domain length
    type: int


steps:
  filter_pfam_structures:
    label: Filter all structural instances for given Pfam families
    run: 'separate_pfam.cwl'
    in:
      track_fams: family_idsfile
      min_dom_len: min_dom_size
    out: [pfam_structs, pfam_obs, splitted_pfam_sep]
     
  resmapping_pfam_structs:
    run:
      class: Workflow
      label: Mapping of residue numbering from UniProt to PDB
      inputs:
        flt_files: 
          type: File
          format: edam:format_3752
        sifts: 
          type: Directory

      steps:
        resmapping_for_Pfam_UP2PDB:
          run: 'resmapping_pfam2pdb.cwl'
          in:
            pfam_sep: flt_files
            sifts_dir: sifts
          out: [ pfam_resmapped, pfam_lost ]

      outputs:
        resmapped_pfam:
          type: File
          format: edam:format_3752
          outputSource: resmapping_for_Pfam_UP2PDB/pfam_resmapped

        lost_insta_pfam:
          type: File
          format: edam:format_2330
          outputSource: resmapping_for_Pfam_UP2PDB/pfam_lost

    scatter: flt_files
    in:
      flt_files: filter_pfam_structures/splitted_pfam_sep
      sifts: siftsdir
    out: [ resmapped_pfam, lost_insta_pfam ]
        
  add_domain_positions:
    label: Add domain positions to residue-mapped instances
    run: 'add_domain_num.cwl'
    in: 
      resmapped_files: resmapping_pfam_structs/resmapped_pfam
      dom_posi_file: resmapped_file
    out: [ resmapped_domains ]

  collect_lost_instances:
    label: Merge lost instances together
    run: 'gather_lost_resmap.cwl'
    in:
      lost_instance: resmapping_pfam_structs/lost_insta_pfam
      obs_insta: filter_pfam_structures/pfam_obs
      outfile: lost_merged
    out: [ lost_domain_list ]


outputs: 
  pfam_domain_posi_file: 
    type: File
    label: All residue-mapped domain instances with domain labels
    format: edam:format_3752
    outputSource: add_domain_positions/resmapped_domains

  pfam_total_lost_structures:
    type: File
    label: Obsolete and inconsistent domain structural instances together 
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
