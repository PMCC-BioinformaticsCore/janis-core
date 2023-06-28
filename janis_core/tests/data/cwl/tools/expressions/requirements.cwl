

#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Tool for testing javascript expression translation / handling

requirements:
  DockerRequirement:
      dockerPull: ubuntu:latest
  InlineJavascriptRequirement: {}
  
  InitialWorkDirRequirement:
    listing: 
    - |
      ${
          return [{"class": "Directory",
                  "basename": "subdir",
                  "listing": [ inputs.example ]
                  }]}
    - entryname: $(1 + 2)
      entry: $("workdir")

  EnvVarRequirement:
    envDef:
    - envName: test1
      envValue: $(9 + 10)

inputs: []
outputs: []

