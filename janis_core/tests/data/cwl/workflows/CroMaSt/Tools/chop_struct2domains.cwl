#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2

label: Compute and write average structure for given set of structures
doc: |
  The tool chops the PDB files into domains according to given residue numebring and average structure from all chopped structures.
  The input file is a json file where key is family id with and value is list of domain structures with residue numbering.

requirements:
  InlineJavascriptRequirement: {}
  EnvVarRequirement:
    envDef:
      KPAX_RESULTS: $(inputs.kpax_result)
  InitialWorkDirRequirement:
    listing:
    - '$({class: "Directory", basename: inputs.pdb_dir.location, listing: []})'
    - |
      ${ 
        if (typeof inputs.split_dir === 'string') {
          return [
            {"class": "Directory", basename: inputs.split_dir, listing: []}]; } 
        else { return [ inputs.split_dir] ; } 
       }
    - entryname: chop_struct2domains.py
      entry:
        $include: Python/chop_struct2domains.py
      writable: false

hints:   
  SoftwareRequirement:
    packages:
      kpax:
        specs: [ "http://kpax.loria.fr/" ]      #link for the resgistration page of this tool is recommended (bio.tools?)
        version: [ "5.1.3.x64" ]


inputs:
  struct_insta:
    label: List of StIs to chop into domains
    type: File
    format: edam:format_3464
    inputBinding:
      position: 1
      prefix: -f

  pdb_dir:
    type: Directory?
    label: The directory containing all the PDB files  
    default: 
      class: Directory
      location: '../PDB_files/'
      basename: PDB_files
    inputBinding:
      position: 2
      prefix: -p

  split_dir:
    type: [ Directory, string, "null"]
    label: The direcory name for storing splitted PDB structures 
    default: split_PDB
    inputBinding:
      position: 3
      prefix: -s

  kpax_result:
    type: string?
    default: KPAX_RESULTS
    inputBinding:
      position: 4
      prefix: -k


outputs:
  split_structs_dir:
    label: Directory with all chopped domains
    type: Directory
    outputBinding:
      glob: $(inputs.split_dir)

  
  family_name:
    label: Instance representing to all chopped domains
    type: stdout


baseCommand:
  - python3
  - chop_struct2domains.py


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
