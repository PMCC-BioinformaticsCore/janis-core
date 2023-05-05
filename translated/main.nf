nextflow.enable.dsl=2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { UNICYCLER } from './modules/unicycler'
include { MULTIQC } from './modules/multiqc'
include { QUAST } from './modules/quast'
include { PROKKA } from './modules/prokka'


// data which will be passed as channels
ch_in_forward_reads  = Channel.fromPath( params.in_forward_reads )
ch_in_long_reads     = Channel.fromPath( params.in_long_reads )
ch_in_reverse_reads  = Channel.fromPath( params.in_reverse_reads )


workflow {

    FASTQC1(
        ch_in_forward_reads
    )

    FASTQC2(
        ch_in_reverse_reads
    )

    UNICYCLER(
        ch_in_forward_reads,
        ch_in_reverse_reads,
        ch_in_long_reads
    )

    MULTIQC(
        FASTQC1.out.outTextFile,
        FASTQC2.out.outTextFile
    )

    QUAST(
        UNICYCLER.out.outAssembly
    )

    PROKKA(
        UNICYCLER.out.outAssembly
    )


}
