#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Add domain position labels to residue-mapped instances
doc: |
  The tool adds domain position labels to each structural instance within the protein in respect with the given list.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.dom_posi_file === 'string') {
          return [
            {"class": "File", "basename": inputs.dom_posi_file, "contents": "", writable: true}]; } 
        else { return [ inputs.dom_posi_file ] ; } 
       }
    - entryname: add_domain_num.py
      entry:
        $include: Python/add_domain_num.py
      writable: false

inputs:
  resmapped_files:
    type: File[]
    format: edam:format_3752
    label: list of files with residue mapped domain structural instances
    inputBinding:
      position: 1
      prefix: -i

  dom_posi_file:
    type: [ File?, string?]
    label: Filename for the output file with domain position labels
    default: resmapped_domains_posi.csv
    inputBinding:
      position: 2
      prefix: -o


outputs:
  resmapped_domains:
    type: File
    label: All residue-mapped domain instances with domain labels
    format: edam:format_3752
    outputBinding: 
      glob: ${ if (typeof inputs.dom_posi_file === 'string') {return inputs.dom_posi_file} else {return [ inputs.dom_posi_file.basename ]}}


baseCommand:
- python3
- add_domain_num.py


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
