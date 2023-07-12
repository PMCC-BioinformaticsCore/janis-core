version development

task Gatk4CalculateContamination {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    File? contaminationTable
    File? statsFile
    File? readOrientationModel
    File pileupTable
    String? segmentationFileOut
    String? contaminationFileOut
  }
  command <<<
    set -e
    gatk CalculateContamination \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if defined(contaminationTable) then ("--contamination-table '" + contaminationTable + "'") else ""} \
      ~{if defined(statsFile) then ("--stats '" + statsFile + "'") else ""} \
      ~{if defined(readOrientationModel) then ("--orientation-bias-artifact-priors '" + readOrientationModel + "'") else ""} \
      -I '~{pileupTable}' \
      --tumor-segmentation '~{select_first([segmentationFileOut, "~{basename(pileupTable)}.mutect2_segments"])}' \
      -O '~{select_first([contaminationFileOut, "~{basename(pileupTable)}.mutect2_contamination"])}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.4.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File contOut = select_first([contaminationFileOut, "~{basename(pileupTable)}.mutect2_contamination"])
    File segOut = select_first([segmentationFileOut, "~{basename(pileupTable)}.mutect2_segments"])
  }
}