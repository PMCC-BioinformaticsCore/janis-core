version development

task Gatk4FilterMutectCalls {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    File? contaminationTable
    File? segmentationFile
    File? statsFile
    File? readOrientationModel
    File vcf
    File vcf_tbi
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputFilename
  }
  command <<<
    set -e
    gatk FilterMutectCalls \
      --java-options '-Xmx~{((select_first([runtime_memory, 16, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if defined(contaminationTable) then ("--contamination-table '" + contaminationTable + "'") else ""} \
      ~{if defined(segmentationFile) then ("--tumor-segmentation '" + segmentationFile + "'") else ""} \
      ~{if defined(statsFile) then ("--stats '" + statsFile + "'") else ""} \
      ~{if defined(readOrientationModel) then ("--orientation-bias-artifact-priors '" + readOrientationModel + "'") else ""} \
      -V '~{vcf}' \
      -R '~{reference}' \
      -O '~{select_first([outputFilename, "generated.vcf.gz"])}'
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
    File out = select_first([outputFilename, "generated.vcf.gz"])
    File out_tbi = select_first([outputFilename, "generated.vcf.gz"]) + ".tbi"
  }
}