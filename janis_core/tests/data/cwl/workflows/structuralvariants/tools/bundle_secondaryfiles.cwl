cwlVersion: v1.0
class: CommandLineTool
id: bundle_secondaryfiles
label: bundle_secondaryfiles

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.primary_file)
      - $(inputs.secondary_files)

baseCommand: [echo]

inputs:
  primary_file:
    type: File
  secondary_files:
    type: 'File[]'

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.primary_file.basename)
    secondaryFiles: ${var arr = []; for (var i = 0; i < inputs.secondary_files.length; i++) { if (inputs.secondary_files[i]) { arr.push(inputs.secondary_files[i].basename) } }; return arr}
