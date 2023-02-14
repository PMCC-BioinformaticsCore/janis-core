cwlVersion: v1.2

class: ExpressionTool

doc: |
  Expression to filter files (by name) in a directory using a regular expression.


requirements:
  InlineJavascriptRequirement: {}
hints:
  LoadListingRequirement:
    loadListing: shallow_listing

inputs:
  folder:
    label: Input folder
    doc: Folder with only files
    type: Directory
  regex:
    label: Regex (JS)
    doc: |
      JavaScript regular expression to be used on the filenames
      MetaBAT2 example: "bin\.[0-9]+\.fa"
    type: string
  output_folder_name:
    type: string
    label: Output folder name
    doc: Output folder name

expression: |
  ${
    var regex = new RegExp(inputs.regex);
    var array = [];
    for (var i = 0; i < inputs.folder.listing.length; i++) {
      if (regex.test(inputs.folder.listing[i].location)){
        array = array.concat(inputs.folder.listing[i]);
      }
    }
    var r = {
       'output_folder':
         { "class": "Directory",
           "basename": inputs.output_folder_name,
           "listing": array
         }
       };
     return r;
  }

outputs:
  output_folder:
    type: Directory

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
