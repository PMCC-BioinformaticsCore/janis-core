
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest

baseCommand: echo

inputs:
 reads: 
  type: File
  inputBinding:
    position: 1
 
 reference: 
  type: File
  inputBinding:
    position: 2
 
 metadata: 
  type: File
  inputBinding:
    position: 3
 
outputs:
  aligned_reads:
    type: stdout