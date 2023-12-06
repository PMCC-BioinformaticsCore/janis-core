
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: ubuntu:latest

baseCommand: echo

inputs:
 the_num_words: 
  type: string
  inputBinding:
    position: 1
 
 the_text: 
  type: string
  inputBinding:
    position: 2

stdout: output.txt

outputs:
  out:
    type: string
    outputBinding:
      glob: output.txt
      loadContents: true
      outputEval: $(self[0].contents.split().length)