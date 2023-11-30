nextflow.enable.dsl=2

include { NGTAX_DEMULTIPLEXING as NGTAX } from './modules/ngtax_demultiplexing'


// data which will be passed as variables
forward_reads  = file( params.forward_reads )
mapping_file   = file( params.mapping_file )
reverse_reads  = file( params.reverse_reads )


workflow {

    // ERROR: PARSING FALLBACKS
    // this is a test
    NGTAX(
        forward_reads,          // forward_reads
        mapping_file,           // mapping_file
        reverse_reads,          // reverse_reads
        params.forward_primer,  // forward_primer
        params.reverse_primer   // reverse_primer
    )


}
