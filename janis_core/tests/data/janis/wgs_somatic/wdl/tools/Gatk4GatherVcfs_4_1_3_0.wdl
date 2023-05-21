version development

task Gatk4GatherVcfs {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    Array[File] vcfs
    String? outputFilename
    Array[File]? argumentsFile
    Int? compressionLevel
    Boolean? createIndex
    Boolean? createMd5File
    File? ga4ghClientSecrets
    Int? maxRecordsInRam
    Boolean? quiet
    File? referenceSequence
    String? tmpDir
    Boolean? useJdkDeflater
    Boolean? useJdkInflater
    String? validationStringency
    Boolean? verbosity
  }
  command <<<
    set -e
    gatk GatherVcfs \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if length(vcfs) > 0 then "--INPUT '" + sep("' --INPUT '", vcfs) + "'" else ""} \
      --OUTPUT '~{select_first([outputFilename, "generated.gathered.vcf"])}' \
      ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
      ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
      ~{if (defined(createIndex) && select_first([createIndex])) then "--CREATE_INDEX" else ""} \
      ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
      ~{if defined(ga4ghClientSecrets) then ("--GA4GH_CLIENT_SECRETS '" + ga4ghClientSecrets + "'") else ""} \
      ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
      ~{if defined(referenceSequence) then ("--REFERENCE_SEQUENCE '" + referenceSequence + "'") else ""} \
      ~{if defined(select_first([tmpDir, "/tmp"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp"]) + "'") else ""} \
      ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--USE_JDK_DEFLATER" else ""} \
      ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--USE_JDK_INFLATER" else ""} \
      ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
      ~{if (defined(verbosity) && select_first([verbosity])) then "--VERBOSITY" else ""}
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
    File out = select_first([outputFilename, "generated.gathered.vcf"])
  }
}