
cwlVersion: v1.0
class: CommandLineTool
 
requirements:
    InlineJavascriptRequirement: {}
    ResourceRequirement:
      coresMin: 8
      ramMin: 3800

baseCommand: [ echo ]

inputs:
    # mandatory
    inFileArr:
        type: File[]
        inputBinding:
            prefix: -i
    inIntArr:
        type: int[]
        inputBinding:
            position: 8
            prefix: "-s"
        doc: Linking distance for stitching
    inStrArr:
        type: string[]
        inputBinding:
            position: 6
        doc: filename to rename to
    
    # optional
    inFileArrOpt:
        type: File[]?
        format:
            - edam:format_1930 # FASTA
            - edam:format_1929 # FASTQ
        inputBinding:
            prefix: -i
    inIntArrOpt:
        type: int[]?
        inputBinding:
            position: 8
            prefix: "-s"
        doc: Linking distance for stitching
    inStrArrOpt:
        type: string[]?
        inputBinding:
            position: 6
        doc: filename to rename to

outputs: []