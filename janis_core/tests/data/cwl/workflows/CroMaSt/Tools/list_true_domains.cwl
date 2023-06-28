#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Format core domain instance list
doc: |
  Fornat core domain instances list from the common instances list identified at first iteration; 
  Preparing input for average structure computation

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.outfile === 'string') {
          return [
            {"class": "File", "basename": inputs.outfile, "contents": "", writable: true}]; } 
        else { return [ inputs.outfile] ; } 
       }
    - entryname: list_true_domains.py
      entry: |-
        #!/usr/bin/env python3
        import json
        all_data = json.load(open('$(inputs.infile.path)', 'r'))
        new = [all_data["0"][a]['$(inputs.database)'] for a in all_data["0"]]
        tmp_dict = {"core": new}
        with open('$(inputs.outfile)', "w") as f:
            json.dump(tmp_dict, f, indent=2)
      writable: false


inputs:
  infile:
    label: True domain StIs from first iteration
    type: File
    format: edam:format_3464

  outfile:
    type: [ File, string, "null"]
    label: User-defined filename for core structures from the given db
    default: coreDomains.json

  database:
    type: [ string, "null"]
    label: The database to select core-structure from; either "CATH" or "Pfam"
    default: CATH


outputs:
  coredomains_list:
    label: Core domain StIs
    type: File
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.outfile === 'string') {return inputs.outfile} else {return inputs.outfile.basename}}


baseCommand:
  - python3
  - list_true_domains.py

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
