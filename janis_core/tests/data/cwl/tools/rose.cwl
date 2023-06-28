cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/rose:v0.0.2
inputs:
  5file:
    type: File
    inputBinding:
      position: 993
      prefix: "--5file"
  annotation_file:
    type: File
    inputBinding:
      position: 7
      prefix: "-g"
    doc: TSV genome annotation file
  bam_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-r"
    secondaryFiles: [".bai"]
    doc: Indexed BAM+BAI file to rank enhancer by
  identifier:
    type: File
    inputBinding:
      position: 990
      prefix: "--identifier"
  metadata:
    type: File
    inputBinding:
      position: 991
      prefix: "--metadata"
  my-file:
    type: File
    inputBinding:
      position: 994
      prefix: "--my-file"
  stitch_distance:
    type: int
    inputBinding:
      position: 8
      prefix: "-s"
    doc: Linking distance for stitching

outputs:
  5plot:
    type: File
    outputBinding:
      glob: "5plot.png"
  e:
    type: File
    outputBinding:
      glob: "e.png"
  gateway_super_enhancers_bed:
    type: File
    outputBinding:
      glob: "*Gateway_SuperEnhancers.bed"
  my-cool-output*:
    type: File
    outputBinding:
      glob: "my-cool-output.gif"
  plot_points_pic:
    type: File
    outputBinding:
      glob: "*Plot_points.png"
baseCommand: ['ROSE_main', '-o', './']
doc: Tool runs ROSE to get Super Enhancers regions