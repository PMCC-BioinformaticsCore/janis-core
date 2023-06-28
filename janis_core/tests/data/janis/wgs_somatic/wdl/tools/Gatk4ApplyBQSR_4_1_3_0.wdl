version development

task Gatk4ApplyBQSR {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputFilename
    File? recalFile
    File? intervals
    Array[String]? intervalStrings
    String? tmpDir
  }
  command <<<
    set -e
    cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
    gatk ApplyBQSR \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      -R '~{reference}' \
      -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' \
      ~{if defined(recalFile) then ("--bqsr-recal-file '" + recalFile + "'") else ""} \
      ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
      ~{if (defined(intervalStrings) && length(select_first([intervalStrings])) > 0) then "--intervals '" + sep("' --intervals '", select_first([intervalStrings])) + "'" else ""} \
      -I '~{bam}' \
      ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--tmp-dir '" + select_first([tmpDir, "/tmp/"]) + "'") else ""}
    if [ -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' ).bai; fi
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
    File out = select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])
    File out_bai = select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"]) + ".bai"
  }
}