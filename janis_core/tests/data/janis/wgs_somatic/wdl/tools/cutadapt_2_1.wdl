version development

task cutadapt {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] fastq
    Array[String]? adapter
    String? outputFilename
    String? secondReadFile
    Int? cores
    String? front
    String? anywhere
    Float? errorRate
    Boolean? noIndels
    Int? times
    Int? overlap
    Boolean? matchReadWildcards
    Boolean? noMatchAdapterWildcards
    String? action
    Int? cut
    String? nextseqTrim
    Int? qualityCutoff
    Boolean? qualityBase
    Int? length
    Int? trimN
    Int? lengthTag
    String? stripSuffix
    String? prefix
    String? suffix
    Boolean? zeroCap
    Int? minimumLength
    Int? maximumLength
    Float? maxN
    Boolean? discardTrimmed
    Boolean? discardUntrimmed
    Boolean? discardCasava
    Boolean? quiet
    String? compressionLevel
    String? infoFile
    String? restFile
    String? wildcardFile
    String? tooShortOutput
    String? tooLongOutput
    String? untrimmedOutput
    Array[String]? removeMiddle3Adapter
    String? removeMiddle5Adapter
    String? removeMiddleBothAdapter
    String? removeNBasesFromSecondRead
    String? pairAdapters
    String? pairFilter
    Boolean? interleaved
    String? untrimmedPairedOutput
    String? tooShortPairedOutput
    String? tooLongPairedOutput
  }
  command <<<
    set -e
    cutadapt \
      ~{if (defined(adapter) && length(select_first([adapter])) > 0) then "-a '" + sep("' -a '", select_first([adapter])) + "'" else ""} \
      -o '~{select_first([outputFilename, "generated--R1.fastq.gz"])}' \
      -p '~{select_first([secondReadFile, "generated--R2.fastq.gz"])}' \
      ~{if defined(cores) then ("--cores " + cores) else ''} \
      ~{if defined(front) then ("--front '" + front + "'") else ""} \
      ~{if defined(anywhere) then ("--anywhere '" + anywhere + "'") else ""} \
      ~{if defined(errorRate) then ("--error-rate " + errorRate) else ''} \
      ~{if (defined(noIndels) && select_first([noIndels])) then "--no-indels" else ""} \
      ~{if defined(times) then ("--times " + times) else ''} \
      ~{if defined(overlap) then ("--overlap " + overlap) else ''} \
      ~{if (defined(matchReadWildcards) && select_first([matchReadWildcards])) then "--match-read-wildcards" else ""} \
      ~{if (defined(noMatchAdapterWildcards) && select_first([noMatchAdapterWildcards])) then "--no-match-adapter-wildcards" else ""} \
      ~{if defined(action) then ("--action '" + action + "'") else ""} \
      ~{if defined(cut) then ("--cut " + cut) else ''} \
      ~{if defined(nextseqTrim) then ("--nextseq-trim '" + nextseqTrim + "'") else ""} \
      ~{if defined(qualityCutoff) then ("--quality-cutoff " + qualityCutoff) else ''} \
      ~{if (defined(qualityBase) && select_first([qualityBase])) then "--quality-base" else ""} \
      ~{if defined(length) then ("--length " + length) else ''} \
      ~{if defined(trimN) then ("--trim-n " + trimN) else ''} \
      ~{if defined(lengthTag) then ("--length-tag " + lengthTag) else ''} \
      ~{if defined(stripSuffix) then ("--strip-suffix '" + stripSuffix + "'") else ""} \
      ~{if defined(prefix) then ("--prefix '" + prefix + "'") else ""} \
      ~{if defined(suffix) then ("--suffix '" + suffix + "'") else ""} \
      ~{if (defined(zeroCap) && select_first([zeroCap])) then "--zero-cap" else ""} \
      ~{if defined(minimumLength) then ("--minimum-length " + minimumLength) else ''} \
      ~{if defined(maximumLength) then ("--maximum-length " + maximumLength) else ''} \
      ~{if defined(maxN) then ("--max-n " + maxN) else ''} \
      ~{if (defined(discardTrimmed) && select_first([discardTrimmed])) then "--discard-trimmed" else ""} \
      ~{if (defined(discardUntrimmed) && select_first([discardUntrimmed])) then "--discard-untrimmed" else ""} \
      ~{if (defined(discardCasava) && select_first([discardCasava])) then "--discard-casava" else ""} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      ~{if defined(compressionLevel) then ("-Z '" + compressionLevel + "'") else ""} \
      ~{if defined(infoFile) then ("--info-file '" + infoFile + "'") else ""} \
      ~{if defined(restFile) then ("--rest-file '" + restFile + "'") else ""} \
      ~{if defined(wildcardFile) then ("--wildcard-file '" + wildcardFile + "'") else ""} \
      ~{if defined(tooShortOutput) then ("--too-short-output '" + tooShortOutput + "'") else ""} \
      ~{if defined(tooLongOutput) then ("--too-long-output '" + tooLongOutput + "'") else ""} \
      ~{if defined(untrimmedOutput) then ("--untrimmed-output '" + untrimmedOutput + "'") else ""} \
      ~{if (defined(removeMiddle3Adapter) && length(select_first([removeMiddle3Adapter])) > 0) then "-A '" + sep("' -A '", select_first([removeMiddle3Adapter])) + "'" else ""} \
      ~{if defined(removeMiddle5Adapter) then ("-G '" + removeMiddle5Adapter + "'") else ""} \
      ~{if defined(removeMiddleBothAdapter) then ("-B '" + removeMiddleBothAdapter + "'") else ""} \
      ~{if defined(removeNBasesFromSecondRead) then ("-U '" + removeNBasesFromSecondRead + "'") else ""} \
      ~{if defined(pairAdapters) then ("--pair-adapters '" + pairAdapters + "'") else ""} \
      ~{if defined(pairFilter) then ("--pair-filter '" + pairFilter + "'") else ""} \
      ~{if (defined(interleaved) && select_first([interleaved])) then "--interleaved" else ""} \
      ~{if defined(untrimmedPairedOutput) then ("--untrimmed-paired-output '" + untrimmedPairedOutput + "'") else ""} \
      ~{if defined(tooShortPairedOutput) then ("--too-short-paired-output '" + tooShortPairedOutput + "'") else ""} \
      ~{if defined(tooLongPairedOutput) then ("--too-long-paired-output '" + tooLongPairedOutput + "'") else ""} \
      ~{if length(fastq) > 0 then "'" + sep("' '", fastq) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 5, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/cutadapt:2.1--py37h14c3975_0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4, 4])}G"
    preemptible: 2
  }
  output {
    Array[File] out = glob("*.fastq.gz")
  }
}