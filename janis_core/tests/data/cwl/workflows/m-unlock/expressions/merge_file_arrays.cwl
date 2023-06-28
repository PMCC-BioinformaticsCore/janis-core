
cwlVersion: v1.2
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement

label: Merge file arrays
doc: |
  Merges arrays of files in an array to a array of files

inputs:
  input:
    type:
      type: array
      items:
        type: array
        items: File

outputs:
  output:
    type: File[]

expression: |
  ${
    var output = [];
    for (var i = 0; i < inputs.input.length; i++) {
      var readgroup_array = inputs.input[i];
      for (var j = 0; j < readgroup_array.length; j++) {
        var readgroup = readgroup_array[j];
        output.push(readgroup);
      }
    }
    return {'output': output}
  }