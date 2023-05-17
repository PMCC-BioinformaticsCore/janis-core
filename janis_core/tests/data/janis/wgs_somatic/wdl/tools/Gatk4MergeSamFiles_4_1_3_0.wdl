version development

task Gatk4MergeSamFiles {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    Array[File] bams
    Array[File] bams_bai
    String? sampleName
    String? outputFilename
    Array[File]? argumentsFile
    Boolean? assumeSorted
    Array[String]? comment
    Boolean? mergeSequenceDictionaries
    String? sortOrder
    Boolean? useThreading
    Int? compressionLevel
    Boolean? createIndex
    Boolean? createMd5File
    Int? maxRecordsInRam
    Boolean? quiet
    File? reference
    File? reference_fai
    File? reference_amb
    File? reference_ann
    File? reference_bwt
    File? reference_pac
    File? reference_sa
    File? reference_dict
    String? tmpDir
    Boolean? useJdkDeflater
    Boolean? useJdkInflater
    String? validationStringency
    String? verbosity
  }
  command <<<
    set -e
    gatk MergeSamFiles \
      --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if (defined(assumeSorted) && select_first([assumeSorted])) then "-AS" else ""} \
      ~{if (defined(comment) && length(select_first([comment])) > 0) then "-CO '" + sep("' '", select_first([comment])) + "'" else ""} \
      ~{if (defined(mergeSequenceDictionaries) && select_first([mergeSequenceDictionaries])) then "-MSD" else ""} \
      ~{if (defined(useThreading) && select_first([useThreading])) then "--USE_THREADING" else ""} \
      ~{if length(bams) > 0 then "-I '" + sep("' -I '", bams) + "'" else ""} \
      -O '~{select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"])}' \
      ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
      ~{if defined(sortOrder) then ("-SO '" + sortOrder + "'") else ""} \
      ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
      ~{if (defined(createIndex) && select_first([createIndex])) then "--CREATE_INDEX" else ""} \
      ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
      ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
      ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
      ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
      ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
      ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
      ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
      ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
    if [ -f $(echo '~{select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"])}' ).bai; fi
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
    File out = select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"])
    File out_bai = select_first([outputFilename, "~{if defined(sampleName) then sampleName else "generated"}.merged.bam"]) + ".bai"
  }
}