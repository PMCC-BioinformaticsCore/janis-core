cwlVersion: v1.2

class: ExpressionTool

doc: |
  Transforms the input folder to an array of files

requirements:
  InlineJavascriptRequirement: {}
hints:
  LoadListingRequirement:
    loadListing: shallow_listing

inputs:
  folder: 
    type: Directory

expression: |
  ${
    var files = [];
    for (var i = 0; i < inputs.folder.listing.length; i++) {
      files.push(inputs.folder.listing[i]);
    }
    return {"files": files};
  }  

outputs:
  files:
    type: File[]

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
s:dateCreated: "2022-10-00"
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
