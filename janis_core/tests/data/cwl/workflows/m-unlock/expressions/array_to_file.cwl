cwlVersion: v1.2
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement

label: Convert an array of 1 file to a file object
doc: |
  Converts the array and returns the first file in the array. 
  Should only be used when 1 file is in the array.

inputs:
  files: 
    type: File[]

outputs:
  file:
    type: File

expression: |
  ${
    var first_file = inputs.files[0];
    return {'file': first_file}
  }