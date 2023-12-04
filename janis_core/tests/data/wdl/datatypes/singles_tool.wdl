version 1.0

task SingleDatatypes {
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

    command <<<
        set -e
        echo \
        
        ~{inFile} \
        ~{inInt} \
        ~{inFloat} \
        ~{inStr} \
        ~{true="--flag1" false="" inBool} \
        
        ~{inFileOpt} \
        ~{inIntOpt} \
        ~{inFloatOpt} \
        ~{inStrOpt} \
        ~{true="--flag2" false="" inBoolOpt} \
        > stdout.txt
    >>>

    output {}

}


