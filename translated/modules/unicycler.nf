nextflow.enable.dsl=2

process UNICYCLER {
    
    container "quay.io/biocontainers/unicycler:0.4.8--py39h98c8e45_5"
    publishDir "${params.outdir}/unicycler"

    input:
    path option11
    path option12
    path option_l, stageAs: 'option_l/*'

    output:
    path "assembly.fasta", emit: outAssembly
    path "assembly.gfa", emit: outAssemblyGraph

    script:
    def option_l = option_l.simpleName != params.NULL ? "-l ${option_l}" : ""
    """
    unicycler \
    -1 ${option11} \
    -2 ${option12} \
    ${option_l} \
    """

}
