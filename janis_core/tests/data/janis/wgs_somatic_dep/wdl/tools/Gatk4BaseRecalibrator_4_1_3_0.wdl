version development

task Gatk4BaseRecalibrator {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    String? tmpDir
    File bam
    File bam_bai
    Array[File] knownSites
    Array[File] knownSites_tbi
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputFilename
    File? intervals
    Array[String]? intervalStrings
  }
  command <<<
    set -e
    cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
    gatk BaseRecalibrator \
      --java-options '-Xmx~{((select_first([runtime_memory, 16, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--tmp-dir '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
      ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
      ~{if (defined(intervalStrings) && length(select_first([intervalStrings])) > 0) then "--intervals '" + sep("' --intervals '", select_first([intervalStrings])) + "'" else ""} \
      -R '~{reference}' \
      -I '~{bam}' \
      -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.table"])}' \
      ~{if length(knownSites) > 0 then "--known-sites '" + sep("' --known-sites '", knownSites) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.3.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 16, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "~{basename(bam, ".bam")}.table"])
  }
}