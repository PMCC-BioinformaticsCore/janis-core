cwlVersion: v1.0
class: CommandLineTool
id: bwa_mem_unpaired
label: bwa_mem_unpaired

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reads)
      - $(inputs.reference_genome)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_5'

baseCommand: [bwa, mem]

inputs:
  reads:
    type: File
    inputBinding:
      position: 3
  reference_genome:
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  threads:
    type: int?
    default: 24
    inputBinding:
      position: 1
      prefix: '-t'
  read_group:
    type: string?
    default: '@RG\tID:SRR709972\tSM:NA06985\tPL:ILLUMINA\tCN:CBRA\tLB:Fragment'
    inputBinding:
      position: 4
      prefix: '-R'

stdout: $(inputs.reads.nameroot.split('.')[0].replace('_', '.R')).sam

outputs:
  output:
    type: stdout
