#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Moves passed files to given directory
doc: |
  The tool copy all the files from their original location to the directory provided by user.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.dir_cp === 'string') {
          return [
            {"class": "Directory", basename: inputs.dir_cp, listing: []}]; } 
        else { return [ inputs.dir_cp] ; } 
       }
    - entryname: script.py
      entry: |-
        #!/usr/bin/env python3
        import json, os
        for i in $(inputs.avg_unp_dom):
            if not isinstance(i, dict):
                continue
            cmd = 'cp {0} {1}'.format(i['path'], '$(inputs.dir_cp)')
            os.system(cmd)

      writable: false

inputs:
  avg_unp_dom:
    label: Files to move into new direcory
    type: File[]
    format: edam:format_1476

  dir_cp:
    type: [ Directory, string, "null"]
    label: The direcory name for storing avg_split PDB structures 
    default: avg_split_PDB

outputs:
  dir_unp_dom:
    type: Directory
    outputBinding:
      glob: $(inputs.dir_cp)
  
baseCommand:
- python
- script.py

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
