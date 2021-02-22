version development

task bwamem {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File reference
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    Array[File] reads
    Array[File]? mates
    String? outputFilename
    Int? threads
    Int? minimumSeedLength
    Int? bandwidth
    Int? offDiagonalXDropoff
    Float? reseedTrigger
    Int? occurenceDiscard
    Boolean? performSW
    Int? matchingScore
    Int? mismatchPenalty
    Int? openGapPenalty
    Int? gapExtensionPenalty
    Int? clippingPenalty
    Int? unpairedReadPenalty
    Boolean? assumeInterleavedFirstInput
    String? readGroupHeaderLine
    Int? outputAlignmentThreshold
    Boolean? outputAllElements
    Boolean? appendComments
    Boolean? hardClipping
    Boolean? markShorterSplits
    Int? verboseLevel
  }
  command <<<
    set -e
    bwa mem \
      ~{if defined(threads) then ("-t " + threads) else ''} \
      ~{if defined(minimumSeedLength) then ("-k " + minimumSeedLength) else ''} \
      ~{if defined(bandwidth) then ("-w " + bandwidth) else ''} \
      ~{if defined(offDiagonalXDropoff) then ("-d " + offDiagonalXDropoff) else ''} \
      ~{if defined(reseedTrigger) then ("-r " + reseedTrigger) else ''} \
      ~{if defined(occurenceDiscard) then ("-c " + occurenceDiscard) else ''} \
      ~{if (defined(performSW) && select_first([performSW])) then "-P" else ""} \
      ~{if defined(matchingScore) then ("-A " + matchingScore) else ''} \
      ~{if defined(mismatchPenalty) then ("-B " + mismatchPenalty) else ''} \
      ~{if defined(openGapPenalty) then ("-O " + openGapPenalty) else ''} \
      ~{if defined(gapExtensionPenalty) then ("-E " + gapExtensionPenalty) else ''} \
      ~{if defined(clippingPenalty) then ("-L " + clippingPenalty) else ''} \
      ~{if defined(unpairedReadPenalty) then ("-U " + unpairedReadPenalty) else ''} \
      ~{if (defined(assumeInterleavedFirstInput) && select_first([assumeInterleavedFirstInput])) then "-p" else ""} \
      ~{if defined(readGroupHeaderLine) then ("-R '" + readGroupHeaderLine + "'") else ""} \
      ~{if defined(outputAlignmentThreshold) then ("-T " + outputAlignmentThreshold) else ''} \
      ~{if (defined(outputAllElements) && select_first([outputAllElements])) then "-a" else ""} \
      ~{if (defined(appendComments) && select_first([appendComments])) then "-C" else ""} \
      ~{if (defined(hardClipping) && select_first([hardClipping])) then "-H" else ""} \
      ~{if (defined(markShorterSplits) && select_first([markShorterSplits])) then "-M" else ""} \
      ~{if defined(verboseLevel) then ("-v " + verboseLevel) else ''} \
      '~{reference}' \
      ~{if length(reads) > 0 then "'" + sep("' '", reads) + "'" else ""} \
      ~{if (defined(mates) && length(select_first([mates])) > 0) then "'" + sep("' '", select_first([mates])) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biocontainers/bwa@sha256:f7b89eccac454a6cf63fac848b98816b0b3a6c857e23f228778bc33b3da2ca2e"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}
