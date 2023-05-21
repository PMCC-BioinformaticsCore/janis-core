nextflow.enable.dsl=2

process BWA_INDEX {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/bwa:0.7.17-dev"
    publishDir "${params.outdir}/prepare_reference/bwa_index"

    input:
    path input_fasta, stageAs: 'input_fasta'

    output:
    path "*.64.alt", optional: true, emit: alt
    path "*.64.amb", optional: true, emit: amb
    path "*.64.ann", optional: true, emit: ann
    path "*.64.bwt", optional: true, emit: bwt
    path "*.64.pac", optional: true, emit: pac
    path "*.64.sa", optional: true, emit: sa

    script:
    """
    \
    <js>inputs.input_alt && inputs.input_amb && inputs.input_ann && inputs.input_bwt && inputs.input_pac && inputs.input_sa ? 'echo bwa' : inputs.generate_bwa_indexes ? 'bwa' : 'echo bwa'</js> \
    index -6 -a bwtsw  \
    ${input_fasta} \
    """

}
