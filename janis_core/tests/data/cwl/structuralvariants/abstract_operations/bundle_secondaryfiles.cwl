cwlVersion: v1.2
class: Operation
id: bundle_secondaryfiles
label: bundle_secondaryfiles

inputs:
  primary_file:
    type: File
  secondary_files:
    type: 'File[]'

outputs:
  output:
    type: File
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
