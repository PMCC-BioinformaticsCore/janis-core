nextflow.enable.dsl=2

include { CUTADAPT } from './modules/cutadapt'


// data which will be passed as channels
ch_in_forward  = Channel.fromPath( params.in_forward )
ch_in_reverse  = Channel.fromPath( params.in_reverse )


workflow {

    CUTADAPT(
        ch_in_forward,  // unknown1
        ch_in_reverse   // unknown2
    )


}
