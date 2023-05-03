nextflow.enable.dsl=2

include { CUTADAPT } from './modules/cutadapt'


// data which will be passed as channels
ch_in_forward   = Channel.fromPath( params.in_forward )
ch_in_reverse   = Channel.fromPath( params.in_reverse )
ch_output_file  = Channel.fromPath( params.output_file )


workflow {

    CUTADAPT(
        ch_output_file,
        ch_in_forward,
        ch_in_reverse
    )


}
