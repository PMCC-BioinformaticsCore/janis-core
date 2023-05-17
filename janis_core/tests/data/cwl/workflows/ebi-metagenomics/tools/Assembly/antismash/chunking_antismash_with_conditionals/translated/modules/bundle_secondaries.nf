nextflow.enable.dsl=2

process BUNDLE_SECONDARIES {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/prepare_reference/bundle_secondaries"

    input:
    path primary_file, stageAs: 'primary_file'
    path secondary_files

    output:
    tuple path("inputs.primary_file.basename"), path("inputs.primary${var arr = []; for (i = 0; i < inputs.secondary_files.length; i++) { if (inputs.secondary_files[i]) { arr.push(inputs.secondary_files[i].basename) } }; return arr}_file.basename"), emit: output

    script:
    """
    echo \
    """

}
