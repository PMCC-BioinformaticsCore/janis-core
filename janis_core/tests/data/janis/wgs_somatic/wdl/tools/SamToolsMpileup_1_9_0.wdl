version development

task SamToolsMpileup {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Boolean? illuminaEncoding
    Boolean? countOrphans
    Boolean? noBAQ
    Int? adjustMQ
    Int? maxDepth
    Boolean? redoBAQ
    File? fastaRef
    File? excludeRG
    File? positions
    Int? minBQ
    Int? minMQ
    String? region
    Boolean? ignoreRG
    String? inclFlags
    String? exclFlags
    Boolean? ignoreOverlaps
    Boolean? outputBP
    Boolean? outputMQ
    Boolean? outputQNAME
    Boolean? allPositions
    Boolean? absolutelyAllPositions
    File? reference
    File bam
    File bam_bai
  }
  command <<<
    set -e
    samtools mpileup \
      ~{if (defined(illuminaEncoding) && select_first([illuminaEncoding])) then "--illumina1.3+" else ""} \
      ~{if (defined(countOrphans) && select_first([countOrphans])) then "--count-orphans" else ""} \
      ~{if (defined(noBAQ) && select_first([noBAQ])) then "--no-BAQ" else ""} \
      ~{if defined(adjustMQ) then ("--adjust-MQ " + adjustMQ) else ''} \
      ~{if defined(maxDepth) then ("--max-depth " + maxDepth) else ''} \
      ~{if (defined(redoBAQ) && select_first([redoBAQ])) then "--redo-BAQ" else ""} \
      ~{if defined(fastaRef) then ("--fasta-ref '" + fastaRef + "'") else ""} \
      ~{if defined(excludeRG) then ("--exclude-RG '" + excludeRG + "'") else ""} \
      ~{if defined(positions) then ("--positions '" + positions + "'") else ""} \
      ~{if defined(minBQ) then ("--min-BQ " + minBQ) else ''} \
      ~{if defined(minMQ) then ("--min-MQ " + minMQ) else ''} \
      ~{if defined(region) then ("--region '" + region + "'") else ""} \
      ~{if (defined(ignoreRG) && select_first([ignoreRG])) then "--ignore-RG" else ""} \
      ~{if defined(inclFlags) then ("--incl-flags '" + inclFlags + "'") else ""} \
      ~{if defined(exclFlags) then ("--excl-flags '" + exclFlags + "'") else ""} \
      ~{if (defined(ignoreOverlaps) && select_first([ignoreOverlaps])) then "--ignore-overlaps" else ""} \
      ~{if (defined(outputBP) && select_first([outputBP])) then "--output-BP" else ""} \
      ~{if (defined(outputMQ) && select_first([outputMQ])) then "--output-MQ" else ""} \
      ~{if (defined(outputQNAME) && select_first([outputQNAME])) then "--output-QNAME" else ""} \
      ~{if (defined(allPositions) && select_first([allPositions])) then "-a" else ""} \
      ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
      '~{bam}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/samtools:1.9--h8571acd_11"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}