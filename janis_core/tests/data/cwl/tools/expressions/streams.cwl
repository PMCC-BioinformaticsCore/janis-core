
#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Tool for testing javascript expression translation / handling

requirements:
  DockerRequirement:
      dockerPull: ubuntu:latest
  InlineJavascriptRequirement: {}

baseCommand: cat

stdin: $(inputs.inFile)
stdout: $(inputs.sampleName.nameroot + ".out")
stderr: $(inputs.sampleName.nameroot + ".err")

inputs:
  inFile: 
    type: File
  sampleName: 
    type: string
    inputBinding:
      position: 1 

outputs: []




