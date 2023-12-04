
version 1.0

import "./io_tool.wdl" as tools

workflow main {
    input {
        Int inInt
        String inStr
        Boolean inBool
        File inFile
        File? inFileOpt
        Array[File] inFileArr
        BamBai inSecondary
    }

    call tools.IoTool {
        input:
            inInt = inInt,
            inStr = inStr,
            inBool = inBool,
            inFile = inFile,
            inFileOpt = inFileOpt,
            inFileArr = inFileArr,
            inSecondary = inSecondary
    }

    output {
        File outFile = IoTool.outFile
        File? outFileOpt = IoTool.outFileOpt
        Array[File] outFileArr = IoTool.outFileArr
        BamBai outSecondary = IoTool.outSecondary
        File outStdout = IoTool.outStdout
    }
}

struct BamBai {
    File bam
    File bai
}