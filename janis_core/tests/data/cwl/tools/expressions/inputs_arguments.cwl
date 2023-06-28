

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}

stdout: output.txt
baseCommand: echo

arguments: 
  - valueFrom: $(runtime.outdir)
    prefix: "-o"
  - valueFrom: $(inputs.inFile.basename)
    prefix: "--name"
  - valueFrom: "--noextract"
  - prefix: -A
    valueFrom: $(1+1)
  - prefix: -B
    valueFrom: $("/foo/bar/baz".split('/').slice(-1)[0])
  - prefix: -C
    valueFrom: |
      ${
        var r = [];
        for (var i = 10; i >= 1; i--) {
          r.push(i);
        }
        return r;
      }

inputs:
  fmt:
    type: string

  inBamBai:
    type: File
    secondaryFiles: $(self.basename + self.nameext.replace('m','i'))
    inputBinding:
      position: 4

  inFile:
    type: File
    inputBinding:
      position: 1

  inFile2:
    type: File
    format: $(inputs.fmt)
    inputBinding:
      position: 2

  runtime_cpu:
    type: int?

  threads:
   type: int
   inputBinding:
     position: 3
     prefix: -t
     valueFrom: |-
       $([inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0])


outputs:
  example_out:
    type: stdout
