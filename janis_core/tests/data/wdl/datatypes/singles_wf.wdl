
version 1.0

import "./singles_tool.wdl" as tools

workflow main {
    input {
        File inFile
        Int inInt
        Float inFloat
        String inStr
        Boolean inBool
        
        File? inFileOpt
        Int? inIntOpt
        Float? inFloatOpt
        String? inStrOpt
        Boolean? inBoolOpt
    }

    call tools.SingleDatatypes {
        input:
            inFile = inFile,
            inInt = inInt,
            inFloat = inFloat,
            inStr = inStr,
            inBool = inBool,

            inFileOpt = inFileOpt,
            inIntOpt = inIntOpt,
            inFloatOpt = inFloatOpt,
            inStrOpt = inStrOpt,
            inBoolOpt = inBoolOpt,
    }

    output {}
}
