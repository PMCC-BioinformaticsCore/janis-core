nextflow.enable.dsl=2

process NGTAX {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/ngtax"

    input:
    path forward_reads, stageAs: 'forward_reads'
    path mapping_file, stageAs: 'mapping_file'
    path reverse_reads, stageAs: 'reverse_reads'

    output:
    path "demultiplexed/*", emit: output

    script:
    """
    -demultiplex -output demultiplexed \
    -mapFile ${mapping_file} \
    -for_p ${params.forward_primer} \
    -rev_p ${params.reverse_primer} \
    -fastQ <js>inputs.forward_reads.path),$(inputs.reverse_reads</js> \
    """

}


process NGTAX {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/ngtax"

    output:
    path "{inputs.sample}_NG-Tax_{inputs.for_read_len}.biom", emit: biom
    path "ngtax2.stdout.log", emit: stdout_out
    path "{inputs.sample}_NG-Tax_{inputs.for_read_len}.ttl", emit: turtle

    script:
    def reference_db = params.reference_db ? "-refdb ${params.reference_db}" : ""
    def rev_read_len = params.rev_read_len ? "-rev_read_len ${params.rev_read_len}" : ""
    def reverse_primer = params.reverse_primer ? "-rev_p ${params.reverse_primer}" : ""
    """
    java -jar /NGTax-2.2.9.jar -ngtax -mapFile cwl_mapping_file.txt \
    ${reference_db} \
    -for_p ${params.forward_primer} \
    -for_read_len ${params.for_read_len} \
    ${reverse_primer} \
    ${rev_read_len} \
    -b "${params.sample}_NG-Tax_${params.for_read_len}.biom" \
    -fragment ${params.fragment} \
    -mock3 ${params.mock3} \
    -mock4 ${params.mock4} \
    -t "${params.sample}_NG-Tax_${params.for_read_len}.ttl" \
    > ngtax2.stdout.log \
    """

}
