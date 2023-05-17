nextflow.enable.dsl=2

process READMAPPING_PILON {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/readmapping_pilon"

    input:
    path forward_reads, stageAs: 'forward_reads'
    path reference, stageAs: 'reference'
    path reverse_reads, stageAs: 'reverse_reads'

    output:
    path "{inputs.identifier}_BBMap_covstats.txt", emit: covstats
    path "{inputs.identifier}_BBMap_log.txt", emit: log
    path "{inputs.identifier}_BBMap.sam", emit: sam
    path "{inputs.identifier}_BBMap_stats.txt", emit: stats

    script:
    def reverse_reads = reverse_reads ? "in2=${reverse_reads}" : ""
    def threads = params.threads ? params.threads : 2
    """
    bbmap.sh \
    threads=${threads} \
    fast=t \
    "-Xmx${params.memory}M" \
    printunmappedcount \
    overwrite=true \
    "statsfile=${params.identifier}_BBMap_stats.txt" \
    "covstats=${params.identifier}_BBMap_covstats.txt" \
    "out=${params.identifier}_BBMap.sam" \
    in=${forward_reads} \
    ${reverse_reads} \
    ref=${reference} \
    2> "${params.identifier}_BBMap_log.txt" \
    """

}
