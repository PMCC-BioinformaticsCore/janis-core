cwlVersion: v1.2
class: ExpressionTool
requirements:
  - class: InlineJavascriptRequirement

label: Merge Directory arrays
doc: |
  Merges arrays of directories in an array to a array of directories

inputs:
  input:
    type:
      type: array
      items:
        type: array
        items: Directory

outputs:
  output:
    type: Directory[]

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