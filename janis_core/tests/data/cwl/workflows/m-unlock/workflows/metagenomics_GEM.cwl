#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
 
label: Metagenomic GEM construction from assembly
doc: |
  Workflow for Metagenomics from bins to metabolic model.<br>
  Summary
    - Prodigal gene prediction
    - CarveMe genome scale metabolic model reconstruction
    - MEMOTE for metabolic model testing
    - SMETANA Species METabolic interaction ANAlysis

  Other UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>
  
  **All tool CWL files and other workflows can be found here:**<br>
    Tools: https://gitlab.com/m-unlock/cwl<br>
    Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>

  **How to setup and use an UNLOCK workflow:**<br>
  https://m-unlock.gitlab.io/docs/setup/setup.html<br>

outputs:
  carveme_gems_folder:
    label: CarveMe GEMs folder
    doc: CarveMe metabolic models folder
    type: Directory
    outputSource: carveme_files_to_folder/results
  protein_fasta_folder:
    label: Protein files folder
    doc: Prodigal predicted proteins (compressed) fasta files
    type: Directory
    outputSource: prodigal_files_to_folder/results
  memote_folder:
    label: MEMOTE outputs folder
    doc: MEMOTE outputs folder
    type: Directory
    outputSource: memote_files_to_folder/results

  smetana_output:
    label: SMETANA output
    doc: SMETANA detailed output table
    type: File?
    outputSource: smetana/detailed_output_tsv

  gemstats_out:
    label: GEMstats
    doc: CarveMe GEM statistics
    type: File
    outputSource: gemstats/carveme_GEMstats

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: Identifier used
  bins:
    type: File[]
    doc: Bin/genome fasta files
    label: Genome/bin
  solver:
    type: string
    doc: Solver to be used in MEMOTE and SMETANA (defaul; cplex)
    default: "cplex"
  run_smetana:
    type: boolean?
    label: Run SMETANA
    doc: Run SMETANA (Species METabolic interaction ANAlysis)
    default: false

  gapfill:
    type: string?
    label: Gap fill
    doc: Gap fill model for given media
    # default: "M8"
  mediadb:
    type: File?
    label: Media database
    doc: Media database file

  threads:
    type: int?
    doc: Number of threads to use for computational processes
    label: number of threads
    default: 2

  destination:
    type: string?
    label: Output Destination (prov only)
    doc: Not used in this workflow. Output destination used for cwl-prov reporting only.

steps:
#############################################
#### Prodigal
  prodigal:
    label: prodigal
    doc: prodigal gene/protein prediction
    run: ../prodigal/prodigal.cwl
    scatter: [input_fasta]
    in:
      input_fasta: bins
      single_mode:
          default: true
    out: [predicted_proteins_faa]

  compress_prodigal:
    label: Compress proteins
    doc: Compress prodigal protein files
    run: ../bash/pigz.cwl
    scatter: inputfile
    in:
      inputfile: 
        source: [prodigal/predicted_proteins_faa]
        linkMerge: merge_flattened
      threads: threads
    out: [outfile]
#############################################
#### CarveMe
  carveme:
    label: CarveMe
    doc: Genome-scale metabolic models reconstruction with CarveMe
    run: ../carveme/carveme.cwl
    scatter: [protein_file]
    in:
      protein_file: prodigal/predicted_proteins_faa
      mediadb: mediadb
      gapfill: gapfill
    out: [carveme_gem]

  compress_carveme:
    label: Compress GEM
    doc: Compress CarveMe GEM
    run: ../bash/pigz.cwl
    scatter: inputfile
    in:
      inputfile: 
        source: [carveme/carveme_gem]
        linkMerge: merge_flattened
      threads: threads
    out: [outfile]
#############################################
#### GEM statistics
  gemstats:
    label: GEM stats
    doc: CarveMe GEM statistics
    run: ../carveme/GEMstats.cwl
    in:
      identifier: identifier
      carveme_gems: carveme/carveme_gem
    out: [carveme_GEMstats]
#############################################
#### SMETANA
  smetana:
    label: SMETANA
    doc: Species METabolic interaction ANAlysis
    when: $(inputs.run_smetana)
    run: ../smetana/smetana.cwl
    in:
      run_smetana: run_smetana
      identifier: identifier
      GEM: carveme/carveme_gem
      solver: solver
    out: [detailed_output_tsv]

#############################################
#### MEMOTE REPORT AND RUN
  memote_report_snapshot:
    label: MEMOTE report snapshot
    doc: Take a snapshot of a model's state and generate a report. 
    run: ../memote/memote.cwl
    scatter: [GEM]
    in:
      GEM: carveme/carveme_gem
      report_snapshot:
        default: true
      skip_test_find_metabolites_produced_with_closed_bounds:
        default: true
      skip_test_find_metabolites_consumed_with_closed_bounds:
        default: true
      skip_test_find_metabolites_not_produced_with_open_bounds:
        default: true
      skip_test_find_metabolites_not_consumed_with_open_bounds:
        default: true
    out: [report_html]

  memote_run:
    label: MEMOTE report snapshot
    doc: MEMOTE run analsis 
    run: ../memote/memote.cwl
    scatter: [GEM]
    in:
      GEM: carveme/carveme_gem
      run:
        default: true
      skip_test_find_metabolites_produced_with_closed_bounds:
        default: true
      skip_test_find_metabolites_consumed_with_closed_bounds:
        default: true
      skip_test_find_metabolites_not_produced_with_open_bounds:
        default: true
      skip_test_find_metabolites_not_consumed_with_open_bounds:
        default: true
    out: [run_json]

#############################################
#### Move to folder if not part of a workflow
  carveme_files_to_folder:
    doc: Preparation of workflow output files to a specific output folder
    label: CarveMe GEMs to folder
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [compress_carveme/outfile]
        linkMerge: merge_flattened
      destination: 
        valueFrom: "CarveMe_GEMs"
    out:
      [results]
#############################################
#### Move to folder if not part of a workflow
  prodigal_files_to_folder:
    doc: Preparation of workflow output files to a specific output folder
    label: Prodigal proteins to folder
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [compress_prodigal/outfile]
        linkMerge: merge_flattened
      destination: 
        valueFrom: "Prodigal_proteins"
    out:
      [results]
#############################################
#### Move to folder if not part of a workflow
  memote_files_to_folder:
    doc: Preparation of workflow output files to a specific output folder
    label: MEMOTE output
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [memote_report_snapshot/report_html,memote_run/run_json]
        linkMerge: merge_flattened
      destination: 
        valueFrom: "MEMOTE"
    out:
      [results]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2022-06-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
