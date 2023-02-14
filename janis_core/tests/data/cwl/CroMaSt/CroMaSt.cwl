#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

inputs:
  pfam: 
    label: Pfam family ids
    type: string[]?
  cath:
    label: CATH family ids
    type: string[]?
  iteration:
    label: Iteration number
    type: int
  filename:
    label: Filename to store family ids per iteration
    type: [ File, string]
  true_domain_file:
    label: To store all the true domain StIs
    type: [ File, string]
  siftsDir:
    label: Directory for storing all SIFTS files 
    type: Directory
  paramfile:
    label: Parameter file for current iteration
    type: File
    format: edam:format_3750
  db_for_core:
    label: Database to select to compute core average structure
    type: string
  core_domain_struct:
    label: Core domain structure (.pdb)
    type: [ File, string]
  prev_crossMapped_pfam:
    label: Pfam cross-mapped domain StIs from previous iteration
    type: File
  prev_crossMapped_cath:
    label: CATH cross-mapped domain StIs from previous iteration
    type: File
  unmapped_analysis_file:
    label: Filename with alignment scores for unmapped instances
    type: string
  pdbDir:
    label: The directory for storing all PDB files
    type: Directory
  cath_resmap:
    label: Filename for residue-mapped CATH domain StIs
    type: string
  cath_lost:
    label: Obsolete and inconsistent CATH domain StIs 
    type: string
  pfam_resmap:
    label: Filename for residue-mapped Pfam domain StIs
    type: string
  pfam_lost:
    label: Obsolete and inconsistent Pfam domain StIs 
    type: string
  domain_like:
    label: To store all the domain-like StIs
    type: [ File, string]
  failed_domain:
    label: To store all failed domain StIs
    type: [ File, string]
  min_domain_length:
    label: Threshold for minimum domain length
    type: int
  alignment_score:
    label: Alignment score from Kpax to analyse structures
    type: string
  score_threshold:
    label: Score threshold for given alignment score from Kpax
    type: float
  unmap_pfam_pass:
    label: Filename to store unmapped but structurally well aligned instances from Pfam
    type: string
  unmap_pfam_fail:
    label: Filename to store unmapped and not properly aligned instances from Pfam
    type: string
  unmap_cath_pass:
    label: Filename to store unmapped but structurally well aligned instances from CATH
    type: string
  unmap_cath_fail:
    label: Filename to store unmapped and not properly aligned instances from CATH
    type: string

requirements:
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

steps:
  get_family_ids:
    label: Get domain family ids
    doc: | 
      Get domain family ids from CATH and Pfam databases from parameter file provided by user
    run: Tools/get_family_ids.cwl
    in:
      pfam_ids: pfam
      cath_ids: cath
      iteration_no: iteration
      fam_tracker: filename
    out: [family_ids]

  pfam_domain_instances:
    label: Produce a list of residue-mapped domain StIs from Pfam ids
    doc: |
      Retrieve and process the PDB structures corresponding to the Pfam family ids resulting in a list of residue-mapped
      structural domain instances along with lost structural instances
      (requires Data/pdbmap downloaded from Pfam and uses SIFTS resource for UniProt to PDB residue Mapping) 
    run: Tools/resmapping_pfam_instances_subwf.cwl
    in:
      family_idsfile: get_family_ids/family_ids
      siftsdir: siftsDir
      resmapped_file: pfam_resmap
      lost_merged: pfam_lost
      min_dom_size: min_domain_length
    out: [ pfam_domain_posi_file, pfam_total_lost_structures ]

  cath_domain_instances:
    label: Produce a list of residue-mapped domain StIs from CATH ids
    doc: |
      Retrieve and process the PDB structures corresponding to the CATH superfamily ids resulting in a list of residue-mapped
      structural domain instances along with lost structural instances
      (requires Data/cath_domain_description_file.txt downloaded from CATH and uses SIFTS resource for PDB to UniProt residue Mapping) 
    run: Tools/resmapping_cath_instances_subwf.cwl
    in:
      family_idsfile: get_family_ids/family_ids
      siftsdir: siftsDir
      resmapped_file: cath_resmap
      lost_merged: cath_lost
      min_dom_size: min_domain_length
    out: [ cath_domain_posi_file, cath_total_lost_structures ]

  add_crossmapped_to_resmapped:
    label: Add cross-mapped to residue-mapped domain StIs
    doc: |
      Add crossmapped domain instances from last iteration to current list of residue mapped domain instances.
    run: Tools/add_crossmapped2resmapped.cwl
    in:
      pfam_resmapped: pfam_domain_instances/pfam_domain_posi_file 
      cath_resmapped: cath_domain_instances/cath_domain_posi_file
      pfam_crossmapped: prev_crossMapped_pfam
      cath_crossmapped: prev_crossMapped_cath
      iteration: iteration
    when: $(inputs.iteration > 0)
    out: [pfam_structs, cath_structs]

  compare_instances_CATH_Pfam:
    label: Compare residue-mapped domain StIs
    doc: |
      Find the intersection between residue-mapped domain StIs of Pfam and CATH lists. 
      Allows variable domain boundaries in a certain range +/- 30aa. Produces three files: common domain instances, 
      and unique domain instances to each Pfam and CATH.
    run: Tools/compare_cath_pfam.cwl
    in:
      iteration: iteration
      resmapped_cath: 
        source: [ cath_domain_instances/cath_domain_posi_file, add_crossmapped_to_resmapped/cath_structs]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }        
      resmapped_pfam: 
        source: [ pfam_domain_instances/pfam_domain_posi_file, add_crossmapped_to_resmapped/pfam_structs]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }
      truedomains_file: true_domain_file
      min_dom_len: min_domain_length
    out: [common_domains, pfam_unique, cath_unique]


  crossmapping_Pfam2CATH:
    label: Map unique Pfam domain StIs to CATH db
    doc: |
      Maps the unique domain StIs from Pfam to the whole CATH database 
      (using residue numbering from PDB allowing variable domain boundaries +/-30aa)
    run: Tools/map_unique_struct_pfam2cath.cwl
    in:
      pfam_unq: compare_instances_CATH_Pfam/pfam_unique
      min_dom_len: min_domain_length
    out: [pfam_crossmapped, allcrossmap_pfam, pfam_unmapped]

  crossmapping_CATH2Pfam:
    label: Map unique CATH domain StIs to Pfam db
    doc: |
      Maps the unique domain StIs from CATH to the whole Pfam database 
      (using residue numbering from UniProt allowing variable domain boundaries +/-30aa)
    run: Tools/map_unique_struct_cath2pfam.cwl
    in:
      cath_unq: compare_instances_CATH_Pfam/cath_unique
      min_dom_len: min_domain_length
    out: [cath_crossmapped, allcrossmap_cath, cath_unmapped]

  format_core_list:
    label: Format core domain StIs list
    doc: |
       Fornat core domain instances list from the common instances list identified at first iteration; 
       Preparing input for average structure computation
    run: Tools/list_true_domains.cwl
    in:
      infile: compare_instances_CATH_Pfam/common_domains
      database: db_for_core
    out: [ coredomains_list ]

  chop_and_avg_for_core:
    label: Compute average of average for core domain instances
    doc: |
      Compute average structure for all averaged structures corresponding to core UniProt domain instances.
      First computes average per UniProt domain instance and then average all averaged structures.
    run: Tools/core_avg_subwf.cwl
    in:
      core_list: format_core_list/coredomains_list
      pdb_dir: pdbDir 
      iteration: iteration
    when: $(inputs.iteration < 1)
    out: [ averaged_structs ]

  chop_and_avg_for_CATH2Pfam:
    label: Compute average of average per cross-mapped Pfam
    doc: |
      Compute average structure for all averaged structures corresponding to UniProt domain instances cross-mapped 
      from CATH to a Pfam family.
      First computes average per UniProt domain instance and then average all averaged structures per Pfam family. 
    run: Tools/other_avg_subwf.cwl
    in:
      crossmap_file: crossmapping_CATH2Pfam/cath_crossmapped
      pdb_dir: pdbDir
    out: [ averaged_structs ]

  chop_and_avg_for_Pfam2CATH:
    label: Compute average of average per cross-mapped CATH
    doc: |
      Compute average structure for all averaged structures corresponding to UniProt domain instances cross-mapped 
      from Pfam to a CATH superfamily.
      First computes average per UniProt domain instance and then average all averaged structures per CATH superfamily. 
    run: Tools/other_avg_subwf.cwl
    in:
      crossmap_file: crossmapping_Pfam2CATH/pfam_crossmapped
      pdb_dir: pdbDir
    out: [ averaged_structs ]

  align_avg_structs_pairwise:
    label: Pairwise alignemnt with core average structure
    doc: |
      Align crossmapped averaged structures against core average domain structure pairwise using Kpax
      Outputs a csv file with all the scores from pairwise alignments
    run: Tools/pairwise_aligner.cwl
    in:
      cath_fam_avg: chop_and_avg_for_CATH2Pfam/averaged_structs
      pfam_fam_avg: chop_and_avg_for_Pfam2CATH/averaged_structs
      core_avg: 
        source: [ chop_and_avg_for_core/averaged_structs, core_domain_struct]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }
      iteration: iteration
    out: [alignment_out] 


  check_alignment_scores:
    label: Checks the alignment score for given threshold
    doc: |
      Checks the alignment score for each aligned structure based on the given threshold
      Outputs the structural instances passing and failing the threshold in separate files 
    run: Tools/check_threshold.cwl
    in:
      aln_res_file: align_avg_structs_pairwise/alignment_out
      fam_tracker: get_family_ids/family_ids
      pfam_crossmap: crossmapping_Pfam2CATH/allcrossmap_pfam
      cath_crossmap: crossmapping_CATH2Pfam/allcrossmap_cath
      aln_score: alignment_score
      threshold_val: score_threshold
    out: [pfam_crossmap_passed, pfam_crossmap_failed, cath_crossmap_passed, cath_crossmap_failed]


  unmapped_from_pfam:
    label: Averages and aligns the unampped instances from Pfam
    doc: |
      First computes average per UniProt domain instance and then aligns all the average structures against core average structure.
      Outputs the alignment results along with the structures passing and failing the threshold for given Kpax score.      
    run: Tools/unmapped_unp_avg_align.cwl
    in:
      unmapped_list: crossmapping_Pfam2CATH/pfam_unmapped
      pdb_dir: pdbDir
      core_struct: 
        source: [ chop_and_avg_for_core/averaged_structs, core_domain_struct]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }
      iteration: iteration
      alignment_score: alignment_score
      score_threshold: score_threshold
      passed_name: unmap_pfam_pass
      failed_name: unmap_pfam_fail
    when: $(inputs.unmapped_list.size > 45)
    out: [ unmapped_aligned_results, domain_like_list, failed_domains_list ]

  unmapped_from_cath:
    label: Averages and aligns the unampped instances from CATH
    doc: |
      First computes average per UniProt domain instance and then aligns all the average structures against core average structure.
      Outputs the alignment results along with the structures passing and failing the threshold for given Kpax score.      
    run: Tools/unmapped_unp_avg_align.cwl
    in:
      unmapped_list: crossmapping_CATH2Pfam/cath_unmapped
      pdb_dir: pdbDir
      core_struct: 
        source: [ chop_and_avg_for_core/averaged_structs, core_domain_struct]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }
      iteration: iteration
      alignment_score: alignment_score
      score_threshold: score_threshold
      passed_name: unmap_cath_pass
      failed_name: unmap_cath_fail
    when: $(inputs.unmapped_list.size > 45)
    out: [ unmapped_aligned_results, domain_like_list, failed_domains_list ]

  gather_domain_like:
    label: Collects all domain-like structural instances 
    doc: |
      Collects all domain-like structural instances from Pfam and CATH
      Outputs the list with all domain-like structural instances together.
    run: Tools/merge_df2json.cwl
    in:
      pfam_unmapped: unmapped_from_pfam/domain_like_list
      cath_unmapped: unmapped_from_cath/domain_like_list
      unmapped_out: domain_like
    out: [ unmapped_list ]

  gather_failed_domains:
    label: Collects all failed domain instances 
    doc: |
      Collects all domain instances failed to pass the criteria from both Pfam and CATH
      Outputs the list with all failed domain instances together.
    run: Tools/merge_df2json.cwl
    in:
      pfam_unmapped: unmapped_from_pfam/failed_domains_list
      cath_unmapped: unmapped_from_cath/failed_domains_list
      pfam_crossmapped: check_alignment_scores/pfam_crossmap_failed
      cath_crossmapped: check_alignment_scores/cath_crossmap_failed
      unmapped_out: failed_domain
    out: [ unmapped_list ]
 

  create_new_parameters:
    label: Create parameter file for next iteration
    doc: |
      Create parameter file for next iteration from previous parameter file
      Filter the pairwise alignments to retrieve family ids passing the threshold for a given Kpax score type       
    run: Tools/create_param.cwl
    in:
      in_paramfile: paramfile
      fam_tracker: get_family_ids/family_ids
      true_domains: compare_instances_CATH_Pfam/common_domains
      core_domain_struct: 
        source: [ chop_and_avg_for_core/averaged_structs, core_domain_struct]
        valueFrom: |
          ${ 
            if (inputs.iteration === 0) { 
              return self[0]; 
            }
            else { return self[1]; } 
          }
      crossmap_pfam: check_alignment_scores/pfam_crossmap_passed
      crossmap_cath: check_alignment_scores/cath_crossmap_passed
      domain_like: gather_domain_like/unmapped_list
      failed_domains: gather_failed_domains/unmapped_list
      iteration: iteration
    out: [next_parmfile]

outputs:
  family_ids_x:
    label: Family ids per iteration
    type: File
    format: edam:format_3464
    outputSource: get_family_ids/family_ids

  resmapped_pfam:
    label: All Pfam residue-mapped domain StIs with domain labels
    type: File
    format: edam:format_3752
    outputSource: pfam_domain_instances/pfam_domain_posi_file

  reslost_pfam:
    label: Obsolete and inconsistent domain StIs from Pfam
    type: File
    format: edam:format_3464
    outputSource: pfam_domain_instances/pfam_total_lost_structures

  resmapped_cath:
    label: All CATH residue-mapped domain StIs with domain labels
    type: File
    format: edam:format_3752
    outputSource: cath_domain_instances/cath_domain_posi_file

  reslost_cath:
    label: Obsolete and inconsistent domain StIs from CATH
    type: File
    format: edam:format_3464
    outputSource: cath_domain_instances/cath_total_lost_structures

  true_domains:
    label: True domain StIs per iteration
    type: File
    format: edam:format_3464
    outputSource: compare_instances_CATH_Pfam/common_domains

  core_domains_list:
    label: Core domain StIs
    type: [ File, "null"]
    format: edam:format_3464
    outputSource: format_core_list/coredomains_list

  core_structure:
    label: Core domain structure (.pdb)
    type: [ File, "null"]
    outputSource: chop_and_avg_for_core/averaged_structs

  all_domain_like:
    label: Domain-like StIs
    type: File
    format: edam:format_3464
    outputSource: gather_domain_like/unmapped_list

  all_failed_domains:
    label: Failed domain StIs
    type: File
    format: edam:format_3464
    outputSource: gather_failed_domains/unmapped_list

  crossmapped_pfam_passed:
    label: Cross-mapped families with Pfam domain StIs passing the threshold
    type: File
    format: edam:format_3464
    outputSource: check_alignment_scores/pfam_crossmap_passed
  
  crossmapped_cath_passed:
    label: Cross-mapped families with CATH domain StIs passing the threshold
    type: File
    format: edam:format_3464
    outputSource: check_alignment_scores/cath_crossmap_passed

  crossres_mappedpfam:
    label: Merged cross-mapped and residue-mapped domain StIs from Pfam
    type: File
    format: edam:format_3752
    outputSource: add_crossmapped_to_resmapped/pfam_structs

  crossres_mappedcath:
    label: Merged cross-mapped and residue-mapped domain StIs from CATH
    type: File
    format: edam:format_3752
    outputSource: add_crossmapped_to_resmapped/cath_structs

  unmap_pfam:
    label: All Pfam un-mapped domin StIs 
    type: File
    format: edam:format_3464
    outputSource: crossmapping_Pfam2CATH/pfam_unmapped

  allmap_pfam:
    label: All Pfam domain StIs cross-mapped to CATH family-wise
    type: File
    format: edam:format_3464
    outputSource: crossmapping_Pfam2CATH/allcrossmap_pfam

  unmap_cath:
    label: All un-mapped domin StIs from CATH 
    type: File
    format: edam:format_3464
    outputSource: crossmapping_CATH2Pfam/cath_unmapped

  allmap_cath:
    label: All CATH cross-mapped domin StIs family-wise together
    type: File
    format: edam:format_3464
    outputSource: crossmapping_CATH2Pfam/allcrossmap_cath

  pfam_crossmap_cath_avg:
    label: Average structures per cross-mapped CATH family for Pfam StIs at family level
    type:
      type: array
      items: File
    format: edam:format_1476
    outputSource: chop_and_avg_for_Pfam2CATH/averaged_structs

  cath_crossmap_pfam_avg:
    label: Average structures per cross-mapped Pfam family for CATH StIs at family level
    type:
      type: array
      items: File
    format: edam:format_1476
    outputSource: chop_and_avg_for_CATH2Pfam/averaged_structs

  avg_alignment_result:
    label: Alignment results from Kpax for all cross-mapped families 
    type: File
    format: edam:format_3752
    outputSource: align_avg_structs_pairwise/alignment_out

  next_parmfile:
    label: Parameter file for next iteration of the workflow
    type: File
    format: edam:format_3750
    outputSource: create_new_parameters/next_parmfile

  align_unmap_pfam:
    label: Alignment results for Pfam unmapped instances
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_pfam/unmapped_aligned_results

  unmap_pfam_passed:
    label: Domain-like StIs from Pfam 
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_pfam/domain_like_list

  unmap_pfam_failed:
    label: Failed domain StIs from Pfam 
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_pfam/failed_domains_list

  align_unmap_cath:
    label: Alignment results for CATH unmapped instances
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_cath/unmapped_aligned_results

  unmap_cath_passed:
    label: Domain-like StIs from CATH 
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_cath/domain_like_list

  unmap_cath_failed:
    label: Failed domain StIs from CATH
    type: File
    format: edam:format_3752
    outputSource: unmapped_from_cath/failed_domains_list

  # unique_pfam:
  #   label: Pfam domain StIs that are not in list of CATH domain StIs
  #   type: File
  #   format: edam:format_3752
  #   outputSource: compare_instances_CATH_Pfam/pfam_unique

  # unique_cath:
  #   label: CATH domain StIs that are not in list of Pfam domain StIs
  #   type: File
  #   format: edam:format_3752
  #   outputSource: compare_instances_CATH_Pfam/cath_unique

  # crossmap_pfam:
  #   label: Pfam domin StIs cross-mapped to CATH family-wise
  #   type: File[]
  #   format: edam:format_3464
  #   outputSource: crossmapping_Pfam2CATH/pfam_crossmapped

  # crossmap_cath:
  #   label: CATH domain StIs cross-mapped to Pfam family-wise
  #   type: File[]
  #   format: edam:format_3464
  #   outputSource: crossmapping_CATH2Pfam/cath_crossmapped  


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
