cwlVersion: v1.2

class: ExpressionTool

doc: |
  Transforms the input files to a mentioned directory

requirements:
 - class: InlineJavascriptRequirement

inputs:
  files:
    type: File[]?
  folders:
    type: Directory[]?
  destination:
    type: string

expression: |
  ${
    var array = []
    if (inputs.files != null) {
      array = array.concat(inputs.files)
    }
    if (inputs.folders != null) {
      array = array.concat(inputs.folders)
    }
    var r = {
       'results':
         { "class": "Directory",
           "basename": inputs.destination,
           "listing": array
         } 
       };
     return r; 
   }

outputs:
  results:
    type: Directory


$namespaces:
 s: http://schema.org/

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"