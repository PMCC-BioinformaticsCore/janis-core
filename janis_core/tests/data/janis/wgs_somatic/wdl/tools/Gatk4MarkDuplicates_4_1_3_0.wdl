version development

task Gatk4MarkDuplicates {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] bam
    String? outputFilename
    String? metricsFilename
    Array[String]? javaOptions
    Int? compression_level
    Array[File]? argumentsFile
    String? assumeSortOrder
    String? barcodeTag
    Array[String]? comment
    Boolean? createIndex
    Boolean? createMd5File
    Int? maxRecordsInRam
    Boolean? quiet
    String? tmpDir
    Boolean? useJdkDeflater
    Boolean? useJdkInflater
    String? validationStringency
    String? verbosity
    Int? opticalDuplicatePixelDistance
  }
  command <<<
    set -e
    gatk MarkDuplicates \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if defined(assumeSortOrder) then ("-ASO '" + assumeSortOrder + "'") else ""} \
      ~{if defined(barcodeTag) then ("--BARCODE_TAG '" + barcodeTag + "'") else ""} \
      ~{if (defined(comment) && length(select_first([comment])) > 0) then "-CO '" + sep("' '", select_first([comment])) + "'" else ""} \
      ~{if defined(opticalDuplicatePixelDistance) then ("--OPTICAL_DUPLICATE_PIXEL_DISTANCE " + opticalDuplicatePixelDistance) else ''} \
      ~{if length(bam) > 0 then "-I '" + sep("' '", bam) + "'" else ""} \
      -O '~{select_first([outputFilename, "generated.markduped.bam"])}' \
      -M '~{select_first([metricsFilename, "generated.metrics.txt"])}' \
      ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
      ~{if select_first([createIndex, true]) then "--CREATE_INDEX" else ""} \
      ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
      ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
      ~{if defined(select_first([tmpDir, "tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "tmp/"]) + "'") else ""} \
      ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
      ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
      ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
      ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
    if [ -f $(echo '~{select_first([outputFilename, "generated.markduped.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "generated.markduped.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "generated.markduped.bam"])}' ).bai; fi
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 4, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.3.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.markduped.bam"])
    File out_bai = select_first([outputFilename, "generated.markduped.bam"]) + ".bai"
    File metrics = select_first([metricsFilename, "generated.metrics.txt"])
  }
}