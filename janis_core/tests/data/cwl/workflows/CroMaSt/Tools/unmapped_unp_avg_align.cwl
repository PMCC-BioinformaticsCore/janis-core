#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: Averages and aligns the unampped instances 
doc: |
  First computes average per UniProt domain instance and then aligns all the average structures against core average structure.
  Outputs the alignment results along with the structures passing and failing the threshold for given Kpax score.      

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  unmapped_list: 
    label: List of un-mapped domin StIs 
    type: File
    format: edam:format_3464
  core_struct:
    label: Core average structure
    type: File
    format: edam:format_1476
  iteration: 
    label: Iteration number
    type: int
  pdb_dir:
    label: The directory with (or to store) PDB files
    type: Directory
  alignment_score:
    label: Alignment score from Kpax to analyse structures
    type: string
  score_threshold:
    label:  Score threshold for given alignment score from Kpax
    type: float
  passed_name:
    label: To store StIs passing the threshold 
    type: string
  failed_name:
    label: To store StIs failed to pass the threshold
    type: string


steps:
  per_unp_dom_instance: #per_unp_dom_instance
    run: 'crossmapped_per_unp_dom.cwl'
    in: 
      fam_structs: unmapped_list
    out: [dom_per_fam, family_name]
        
  avg_unp_domains:
    run:
      class: Workflow
      inputs:
        domfiles: File
        pdb_d: Directory
      outputs:
        avg_unp_dom_structs:
          type: File
          format: edam:format_1476
          outputSource: avg_chopped_structs_unp_domains/avg_structs
            
      steps:
        chop_structs:
          run: 'chop_struct2domains.cwl'
          in:
            struct_insta: domfiles
            pdb_dir: pdb_d
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

    scatter: domfiles
    in:
      domfiles: per_unp_dom_instance/dom_per_fam
      pdb_d: pdb_dir
    out: [ avg_unp_dom_structs ]
        
  copy_avg_dom:
    run: 'move_files.cwl'
    in: 
      avg_unp_dom: avg_unp_domains/avg_unp_dom_structs
    out: [ dir_unp_dom ]

  pairwise_align_avg_structs:
    run: 'pairwise_aligner.cwl'
    in:
      query_dir: copy_avg_dom/dir_unp_dom
      core_avg: core_struct
    out: [alignment_out]
  
  check_threshold_step:
    run: 'filter_align_scores.cwl'
    in:
      aln_result: pairwise_align_avg_structs/alignment_out
      all_struct_list: unmapped_list
      aln_score: alignment_score
      threshold_val: score_threshold
      passed_structs: passed_name
      failed_structs: failed_name
    out: [ passed_structs_list, failed_structs_list]


outputs: 
  unmapped_aligned_results: 
    label: Alignment results from Kpax for unmapped instances
    type: File
    format: edam:format_3752
    outputSource: pairwise_align_avg_structs/alignment_out
  
  domain_like_list:
    label: Domain-like StIs
    type: File
    format: edam:format_3752
    outputSource: check_threshold_step/passed_structs_list

  failed_domains_list:
    label: Failed domain StIs
    type: File
    outputSource: check_threshold_step/failed_structs_list
    format: edam:format_3752

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
