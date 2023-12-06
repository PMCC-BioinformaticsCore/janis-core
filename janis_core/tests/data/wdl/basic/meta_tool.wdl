version 1.0

task Rename {
    runtime {
        docker: "biowardrobe2/scidap:v0.0.3"
    }

    input {
        File sourceFile
        String targetFilename
    }

    command <<<
        set -e
        cp \
        ~{sourceFile} \
        ~{targetFilename} \
    >>>

    output {
        File targetFile = targetFilename
    }

    parameter_meta {
        sourceFile: {
            description: "source file to rename"
        },
        targetFilename: {
            description: "filename to rename to"
        }
    }
}
