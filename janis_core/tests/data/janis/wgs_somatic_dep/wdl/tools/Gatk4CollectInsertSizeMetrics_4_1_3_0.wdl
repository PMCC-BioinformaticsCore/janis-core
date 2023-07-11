version development

task Gatk4CollectInsertSizeMetrics {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    File bam
    File bam_bai
    String? outputFilename
    String? outputHistogram
    Array[File]? argumentsFile
    Boolean? assumeSorted
    Float? deviations
    Int? histogramWidth
    Boolean? includeDuplicates
    String? metricAccumulationLevel
    Float? minimumPCT
    Int? stopAfter
    Boolean? version
    Boolean? showHidden
  }
  command <<<
    set -e
    gatk CollectInsertSizeMetrics \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      -O '~{select_first([outputFilename, "generated.metrics.txt"])}' \
      -H '~{select_first([outputHistogram, "generated.histogram.pdf"])}' \
      -I '~{bam}' \
      ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' --arguments_file '", select_first([argumentsFile])) + "'" else ""} \
      ~{if (defined(assumeSorted) && select_first([assumeSorted])) then "--ASSUME_SORTED" else ""} \
      ~{if defined(deviations) then ("--DEVIATIONS " + deviations) else ''} \
      ~{if defined(histogramWidth) then ("--HISTOGRAM_WIDTH " + histogramWidth) else ''} \
      ~{if (defined(includeDuplicates) && select_first([includeDuplicates])) then "--INCLUDE_DUPLICATES" else ""} \
      ~{if defined(metricAccumulationLevel) then ("--METRIC_ACCUMULATION_LEVEL '" + metricAccumulationLevel + "'") else ""} \
      ~{if defined(minimumPCT) then ("--MINIMUM_PCT " + minimumPCT) else ''} \
      ~{if defined(stopAfter) then ("--STOP_AFTER " + stopAfter) else ''} \
      ~{if (defined(version) && select_first([version])) then "--version" else ""} \
      ~{if (defined(showHidden) && select_first([showHidden])) then "--showHidden" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.3.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.metrics.txt"])
    File outHistogram = select_first([outputHistogram, "generated.histogram.pdf"])
  }
}