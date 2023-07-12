#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "test workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    inFile:
        type: File
    inFileOptional:
        type: File?
    inFileArr:
        type: File[]
    inFileArrOptional:
        type: File[]?
    inSecondary:
        type: File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    inString:
        type: string
    inStringOptional:
        type: string?
    inStringOptional2:
        type: string?
    inStringArr:
        type: string[]
    inStringArrOptional:
        type: string[]?
    inInt:
        type: int
    inIntOptional:
        type: int?
    inIntOptional2:
        type: int?
        default: 5

outputs:
    outFile:
        type: File
        outputSource: optional1/out_stdout

steps:
    basic:
        run: basic.cwl
        in:
            inFile: inFile
            inString: inString
            inSecondary: inSecondary
            inFileOptional: inFileOptional
            inIntOptional: inIntOptional2
        out:
            [out_stdout]
    mandatory:
        run: mandatory_input_types.cwl
        in:
            inFile: inFile
            inFileArr: inFileArr
            inSecondary: inSecondary
            inString: inString
            inStringArr: inStringArr
            inInt: inInt
        out:
            [out_stdout]
    optional1:
        run: optional_input_types.cwl
        in:
            inFile: inFileOptional
            inFileArr: inFileArrOptional
        out:
            [out_stdout]
    optional2:
        run: optional_input_types.cwl
        in:
            inString: inStringOptional
            inStringArr: inStringArrOptional
            inInt: inIntOptional
        out:
            [out_stdout]
    sub:
        run: subworkflow.cwl
        in:
            inFile: inFile
            inFileArr: inFileArr
            inSecondary: inSecondary
            inString: inString
            inStringArr: inStringArr
            inInt: inInt
        out:
            [out_stdout]
