#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: MEMOTE

doc: MEMOTE, short for metabolic model testing   

requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "MEMOTE"
      writable: true
    - entryname: script.sh
      entry: |-
          #!/bin/bash
          source /root/miniconda/bin/activate
          conda init bash
          conda activate /unlock/infrastructure/conda/memote/cplex/memote_0.13.0
          memote $@


baseCommand: ["bash", "script.sh"] # see requirements

inputs:
  GEM:
    type: File
    doc: Metabolic model (xml format)
    label: Metabolic model
    inputBinding:
      position: 100
  report_snapshot:
    type: boolean?
    label: Report
    doc: Take a snapshot of a model's state and generate a report.
    inputBinding:
      position: 0
  run:
    type: boolean?
    label: Run
    doc: Run the test suite on a single model and collect results.
    inputBinding:
      position: 0
  
  solver:
    type: string?
    label: solver
    doc: Set the solver to be used. [cplex|glpk|gurobi|glpk_exact]. default; glpk
    inputBinding:
      prefix: --solver
      position: 3

  skip_test_find_metabolites_produced_with_closed_bounds:
    type: boolean?
    doc: Skip test; find metabolites produced with closed bounds
    label: Skip test; find metabolites produced with closed bounds
    inputBinding:
      prefix: "--skip test_find_metabolites_produced_with_closed_bounds"
      position: 11
  skip_test_find_metabolites_consumed_with_closed_bounds:
    type: boolean?
    doc: Skip test; Find metabolites consumed with closed bounds
    label: Skip test; Find metabolites consumed with closed bounds
    inputBinding:
      prefix: "--skip test_find_metabolites_consumed_with_closed_bounds"
      position: 12
  skip_test_find_metabolites_not_produced_with_open_bounds:
    type: boolean?
    doc: Skip test; Find metabolites not produced with open bounds
    label: Skip test; Find metabolites not produced with open bounds
    inputBinding:
      prefix: "--skip test_find_metabolites_not_produced_with_open_bounds"
      position: 13
  skip_test_find_metabolites_not_consumed_with_open_bounds:
    type: boolean?
    doc: Skip test; Find metabolites not consumedwith open_bounds
    label: Skip test; Find metabolites not consumedwith open_bounds
    inputBinding:
      prefix: "--skip test_find_metabolites_not_consumed_with_open_bounds"
      position: 14
  skip_test_find_incorrect_thermodynamic_reversibility:
    type: boolean?
    doc: Skip test; find incorrect thermodynamic reversibility
    label: Skip test; find incorrect thermodynamic reversibility
    inputBinding:
      prefix: "--skip skip_test_find_incorrect_thermodynamic_reversibility"
      position: 15

arguments:
  - |
    ${
      if (inputs.run){
        return "run --filename " + inputs.GEM.basename + "_MEMOTE.json.gz";
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.report_snapshot){
        return "report snapshot --filename " + inputs.GEM.basename + "_MEMOTE.html";
      } else {
        return '';
      }
    }

outputs:
  run_json:
    type: File?
    outputBinding:
      glob: $(inputs.GEM.basename)_MEMOTE.json.gz
  report_html:
    type: File?
    outputBinding:
      glob: $(inputs.GEM.basename)_MEMOTE.html


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
