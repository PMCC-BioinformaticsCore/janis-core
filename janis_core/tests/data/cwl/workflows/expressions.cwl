

#!/usr/bin/env cwl-runner
class: Workflow
cwlVersion: v1.2

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement

inputs:

  bambai_pair1:
    type: File
    secondaryFiles:
    - .bai
    label: "Coordinate sorted BAM+BAI files"
    format: "http://edamontology.org/format_2572"
    doc: "Coordinate sorted BAM file and BAI index file"

  bambai_pair2:
    type: File
    secondaryFiles:
      pattern: $(self.basename + self.nameext.replace('m','i'))
  
  fmt:
    type: string

  png:
    type: File

  prefix_str:
    type: string
  
  suffix_str:
    type: string

outputs:
  out1:
    type: File
    outputSource: rename_png1/target_file
  
  out2:
    type: File
    format: $(inputs.fmt)
    outputSource: rename_png2/target_file

  out3:
    type: File
    secondaryFiles:
      pattern: $(self.basename + self.nameext.replace('m','i'))
    outputSource: bambai_pair2

steps:
  echo1:
    in:
      text: prefix_str
    out: [out]
    run: ../tools/echo.cwl
  
  echo2:
    in:
      text: 
        valueFrom: $(inputs.prefix_str + '_' + inputs.suffix_str)
    out: [out]
    run: ../tools/echo.cwl
        
  rename_png1:
    in:
      source_file: png
      target_filename:
        source: bambai_pair1
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+"_default_s_enhcr.png")
    out: [target_file]
    run: ../tools/rename_png.cwl
  
  rename_png2:
    in:
      source_file: png
      target_filename:
        source: bambai_pair2
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+"_default_s_enhcr.png")
    out: [target_file]
    run: ../tools/rename_png.cwl


