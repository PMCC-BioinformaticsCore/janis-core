#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Changes the format for core structural instances (only 1st iteration) 
doc: |
  The tool reads the given family IDs from parameter file (.yml) and writes it to a separate file according to each iteration.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - |
      ${ 
        if (typeof inputs.outfile === 'string') {
          return [
            {"class": "File", "basename": inputs.outfile, "contents": "{}", writable: true}]; } 
        else { return [ inputs.outfile] ; } 
       }
    - entryname: collect_lost_instances.py
      entry: |-
        #!/usr/bin/env python3
        import sys, json
        lost_whl = []
        total_content = json.load(open(sys.argv[1], 'r'))
        curr_iter = str(len(total_content))
        for i in $(inputs.lost_instance):
            if not isinstance(i, dict):
                continue
            lost_whl.extend(open(i['path']).read().split('\n'))
        lost_whl = [i for i in lost_whl if i]
        obs_whl = open('$(inputs.obs_insta.path)').read().split('\n')
        obs_whl = [i for i in obs_whl if i]

        total_content[curr_iter] = {"Obsolete": [*set(obs_whl)], "Lost_res-map": [*set(lost_whl)]}
        with open(sys.argv[1], 'w') as f:
            json.dump(total_content, f, indent=2, sort_keys=True)
      writable: false


inputs:
  lost_instance:
    type: File[]
    label: The files containing inconsistent entries found during residue-mapping step
    format: edam:format_2330

  obs_insta:
    type: File
    label: The file containing obsolete domain structural instances for given family
    format: edam:format_2330

  outfile:
    type: [ File, string, "null"]
    label: User-defined filename for lost domain structural instances while resmapping
    default: lost_resmap.json
    inputBinding:
      position: 1


outputs:
  lost_domain_list:
    type: File
    label: Obsolete and inconsistent domain structural instances together 
    format: edam:format_3464
    outputBinding:
      glob: ${ if (typeof inputs.outfile === 'string') {return inputs.outfile} else {return inputs.outfile.basename}}


baseCommand:
  - python3
  - collect_lost_instances.py


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
