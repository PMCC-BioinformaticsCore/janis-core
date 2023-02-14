#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Divides a dictionary of structures per family into multiple according to domain instance
doc: |
  The tool takes a dictionary and divides all the structures according to their UniProt id and domain position.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - entryname: crossmapped_per_unp_dom.py
      entry: |-
        #!/usr/bin/env python3
        import json, sys, os
        all_data = json.load(open('$(inputs.fam_structs.path)', 'r'))
        done_dom = {}
        # This is designed only for first key in the dict/json 
        if len(all_data.keys()) < 1:
            sys.exit()
        fam_name = list(all_data.keys())[0]
        for struct in all_data[fam_name]:
            struct = struct.split(',')
            if fam_name + '_' + '_'.join(struct[6:8]) in done_dom:
                done_dom[fam_name + '_' + '_'.join(struct[6:8])].append(','.join(struct))
            else:
                done_dom[fam_name + '_' + '_'.join(struct[6:8])] = [','.join(struct)]

        for domx in done_dom:
            tmp_dom = {domx: done_dom[domx]}
            with open(domx + '.json', 'w') as f:
                json.dump(tmp_dom, f, indent=2)
        print(fam_name)
      writable: false

inputs:
  fam_structs:
    label: Cross-mapped domain StIs per family
    type: File
    format: edam:format_3464
    inputBinding:
      position: 1
      prefix: -f

outputs:
  dom_per_fam:
    label: List of StIs per UniProt domain instance
    type: 
      type: array
      items: File
    format: edam:format_3464
    outputBinding:
      glob: "*.json"

  family_name:
    label: UniProt domain instance corresponding to StIs
    type: stdout

baseCommand:
- python3
- crossmapped_per_unp_dom.py


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
