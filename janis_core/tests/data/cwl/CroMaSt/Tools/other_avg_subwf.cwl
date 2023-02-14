#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: Compute average of average per cross-mapped famil(y)ies
doc: |
  Compute average structure for all averaged structures corresponding to UniProt domain instances cross-mapped 
  from Pfam/CATH to a CATH/Pfam family.
  First computes average per UniProt domain instance and then average all averaged structures per Pfam family. 

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}


inputs:
  crossmap_file: 
    label: CATH cross-mapped domin StIs family-wise 
    type: File[]
    format: edam:format_3464
  pdb_dir:
    label: The directory with (or to store) PDB files
    type: Directory


steps:
  chop_and_avg_from_list:
    run: 
      class: Workflow
      inputs:
        in_file:
          label: CATH cross-mapped domin StIs per family
          type: File
          format: edam:format_3464
        pdb_storage: 
          label: The directory with (or to store) PDB files
          type: Directory

      steps:
        per_dom_instance:
          run: 'crossmapped_per_unp_dom.cwl'
          in: 
            fam_structs: in_file
          out: [dom_per_fam, family_name]
        
        avg_unp_domains:
          run:
            class: Workflow
            inputs:
              dom_structs: File
              pdbDir: Directory
            
            steps:
              chop_structs:
                run: 'chop_struct2domains.cwl'
                in:
                  struct_insta: dom_structs
                  pdb_dir: pdbDir
                out: [split_structs_dir, family_name]

              avg_chopped_structures_unp_domains:
                run: 'align_compute_avg.cwl'
                in:
                  fam_name: chop_structs/family_name
                  split_dir: chop_structs/split_structs_dir
                out: [avg_structs]

            outputs:
              avg_structures:
                label: Average structure per UniProt domain instance
                type: File
                format: edam:format_1476
                outputSource: avg_chopped_structures_unp_domains/avg_structs

          scatter: dom_structs
          in:
            dom_structs: per_dom_instance/dom_per_fam
            pdbDir: pdb_storage
          out: [ avg_structures ]
        
        copy_avg_dom:
          run: 'move_files.cwl'
          in: 
            avg_unp_dom: avg_unp_domains/avg_structures
          out: [ dir_unp_dom ]

        avg_averaged_unp_domains:
          run: 'align_compute_avg.cwl'
          in:
            fam_name: per_dom_instance/family_name
            split_dir: copy_avg_dom/dir_unp_dom
          out: [avg_structs]

      outputs: 
        avg_struct_per_fam:
          label: Avergae structure for the family with passed StIs
          type: File
          format: edam:format_1476
          outputSource: avg_averaged_unp_domains/avg_structs
        
    scatter: in_file
    in:
      in_file: crossmap_file
      pdb_storage: pdb_dir
    out: [ avg_struct_per_fam ]


outputs:
  averaged_structs:
    label: Average structures per passed family with given StIs
    type:
      type: array
      items: File
    format: edam:format_1476
    outputSource: chop_and_avg_from_list/avg_struct_per_fam
  

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

