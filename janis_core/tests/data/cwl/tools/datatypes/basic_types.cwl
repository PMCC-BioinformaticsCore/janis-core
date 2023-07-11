
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
    inFile:
        type: File
        format:
            - edam:format_1930 # FASTA
            - edam:format_1929 # FASTQ
        inputBinding:
            prefix: -i
    inInt:
        type: int
        inputBinding:
            position: 8
            prefix: "-s"
        doc: Linking distance for stitching
    inStr:
        type: string
        inputBinding:
            position: 6
        doc: filename to rename to
    
    # optional
    inFileOpt:
        type: File?
        format:
            - edam:format_1930 # FASTA
            - edam:format_1929 # FASTQ
        inputBinding:
            prefix: -i
    inIntOpt:
        type: int?
        inputBinding:
            position: 8
            prefix: "-s"
        doc: Linking distance for stitching
    inStrOpt:
        type: string?
        inputBinding:
            position: 6
        doc: filename to rename to

outputs: []
    