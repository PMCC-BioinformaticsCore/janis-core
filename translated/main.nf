nextflow.enable.dsl=2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { UNICYCLER } from './modules/unicycler'
include { NANOPLOT } from './modules/nanoplot'
include { QUAST } from './modules/quast'
include { BUSCO } from './modules/busco'


// data which will be passed as channels
ch_in_long      = Channel.fromPath( params.in_long )
ch_in_short_r1  = Channel.fromPath( params.in_short_r1 )
ch_in_short_r2  = Channel.fromPath( params.in_short_r2 )


workflow {

    FASTQC1(
        ch_in_short_r1  // inputFile
    )

    FASTQC2(
        ch_in_short_r2  // inputFile
    )

    UNICYCLER(
        ch_in_short_r1,  // option11
        ch_in_short_r2,  // option12
        ch_in_long       // optionL
    )

    NANOPLOT(
        ch_in_long  // unknown1
    )

    QUAST(
        UNICYCLER.out.outAssembly  // unknown1
    )

    BUSCO(
        UNICYCLER.out.outAssembly  // inFile
    )


}
