
class: CommandLineTool
cwlVersion: v1.0
doc: "Reverse each line using the `rev` command"

hints:
  ResourceRequirement:
    ramMin: 8
  DockerRequirement:
    dockerPull: "debian:stretch-slim"

baseCommand: rev
stdout: output.txt

inputs:
  bam:
    type: File
  
  fastq:
    type: File
  
  infile:
    type: File
    inputBinding: 
      position: 1
    format: edam:format_2330

outputs:
  out1:
    type: File
    outputBinding:
      glob: output.txt
    format: $(inputs.infile.format)

  out2:
    type: File
    outputBinding:
      glob: output.txt
    secondaryFiles: $(inputs.bam.basename + inputs.bam.nameext.replace('m','i'))

  out3:
    type: File
    outputBinding:
      glob: $(inputs.fastq.nameroot.replace(/\b.fastq\b/g, '')).trimmed.fastq
  
  out4:
    type: File
    outputBinding:
      outputEval: $(self[0].contents)

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - EDAM.owl
