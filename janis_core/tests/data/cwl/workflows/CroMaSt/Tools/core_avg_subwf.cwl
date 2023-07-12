#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: Compute average of average for core domain instances
doc: |
  Compute average structure for all averaged structures corresponding to core UniProt domain instances.
  First computes average per UniProt domain instance and then average all averaged structures.

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  core_list:
    label: Core domain StIs 
    type: File
    format: edam:format_3464
  pdb_dir:
    label: The directory with (or to store) PDB files
    type: Directory


steps:
  per_unp_dom_instance: 
    run: 'crossmapped_per_unp_dom.cwl'
    in: 
      fam_structs: core_list
    out: [dom_per_fam, family_name]
        
  avg_unp_domains:
    run:
      class: Workflow
      inputs:
        domfiles: File
        pdb_storage: Directory

      steps:
        chop_structs:
          run: 'chop_struct2domains.cwl'
          in:
            struct_insta: domfiles
            pdb_dir: pdb_storage
          out: [split_structs_dir, family_name]

        avg_chopped_structs_unp_domains:
          label: Compute average structure for each UniProt domain
          doc: |
            Align all the provided structures using Kpax. And compute an average structure from all the aligned structural instances,
            based on aligned residues for at least 50% of the instances. 
          run: 'align_compute_avg.cwl'
          in:
            fam_name: chop_structs/family_name
            split_dir: chop_structs/split_structs_dir
          out: [avg_structs]

      outputs:
        avg_unp_dom_structs:
          label: Computed average structure for passed structures
          type: File
          format: edam:format_1476
          outputSource: avg_chopped_structs_unp_domains/avg_structs

    scatter: domfiles
    in:
      domfiles: per_unp_dom_instance/dom_per_fam
      pdb_storage: pdb_dir
    out: [ avg_unp_dom_structs ]
        
  copy_avg_dom:
    run: 'move_files.cwl'
    in: 
      avg_unp_dom: avg_unp_domains/avg_unp_dom_structs
    out: [ dir_unp_dom ]

  avg_averaged_structures: 
    label: Compute average structure for all provided averaged structures
    doc: |
      Align all the provided averaged structures using Kpax. And compute an average structure from all the aligned averaged structures,
      based on aligned residues for at least 50% of the structures. 
    run: 'align_compute_avg.cwl'
    in:
      fam_name: per_unp_dom_instance/family_name
      split_dir: copy_avg_dom/dir_unp_dom
    out: [avg_structs]


outputs: 
  averaged_structs: 
    label: Core average structure
    type: File
    format: edam:format_1476
    outputSource: avg_averaged_structures/avg_structs
  
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
