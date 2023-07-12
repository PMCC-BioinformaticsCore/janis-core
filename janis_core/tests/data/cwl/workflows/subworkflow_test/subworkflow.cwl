#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "test subworkflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    inFile:
        type: File
    inFileArr:
        type: File[]
    inSecondary:
        type: File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    inString:
        type: string
    inStringArr:
        type: string[]
    inInt:
        type: int

outputs:
    outFile:
        type: File
        outputSource: optional/out_stdout

steps:
    optional:
        run: optional_input_types.cwl
        in:
            inFile: inFile
            inFileArr: inFileArr
        out:
            [out_stdout]
