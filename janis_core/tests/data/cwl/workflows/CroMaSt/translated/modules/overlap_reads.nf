nextflow.enable.dsl=2

process OVERLAP_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/overlap_reads"
    cpus "${params.before_qc.overlap_reads.overlap_reads.cpus}"
    memory "${params.before_qc.overlap_reads.overlap_reads.memory}"

    input:
    path forward_reads, stageAs: 'forward_reads'
    path namefile, stageAs: 'namefile'
    path reverse_reads, stageAs: 'reverse_reads'

    output:
    path "forward_unmerged.fastq.gz", emit: forward_unmerged_reads
    path "*_MERGED*", emit: merged_reads
    path "reverse_unmerged.fastq.gz", emit: reverse_unmerged_reads

    script:
    def forward_reads = forward_reads ? "-f ${forward_reads}" : ""
    def reverse_reads = reverse_reads ? "-r ${reverse_reads}" : ""
    """
    SeqPrep \
    ${forward_reads} \
    ${reverse_reads} \
    -s <js>${ return inputs.namefile.nameroot.split('_')[0] + '_MERGED.fastq.gz' }</js> \
    -1 \
    forward_unmerged.fastq.gz \
    -2 \
    reverse_unmerged.fastq.gz \
    """

}


process OVERLAP_READS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/before_qc/overlap_reads/overlap_reads"
    cpus "${params.before_qc.overlap_reads.overlap_reads.cpus}"
    memory "${params.before_qc.overlap_reads.overlap_reads.memory}"

    input:
    path forward_reads, stageAs: 'forward_reads'
    path namefile, stageAs: 'namefile'
    path reverse_reads, stageAs: 'reverse_reads'

    output:
    path "forward_unmerged.fastq.gz", emit: forward_unmerged_reads
    path "*_MERGED*", emit: merged_reads
    path "reverse_unmerged.fastq.gz", emit: reverse_unmerged_reads

    script:
    def forward_reads = forward_reads ? "-f ${forward_reads}" : ""
    def reverse_reads = reverse_reads ? "-r ${reverse_reads}" : ""
    """
    SeqPrep \
    ${forward_reads} \
    ${reverse_reads} \
    -s <js>${ return inputs.namefile.nameroot.split('_')[0] + '_MERGED.fastq.gz' }</js> \
    -1 \
    forward_unmerged.fastq.gz \
    -2 \
    reverse_unmerged.fastq.gz \
    """

}
