version development

task SamToolsView {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Boolean? cramOutput
    Boolean? compressedBam
    Boolean? uncompressedBam
    Boolean? onlyOutputHeader
    Boolean? countAlignments
    File? writeAlignments
    File? inputTSV
    File? onlyOverlapping
    Boolean? useMultiRegionIterator
    String? outputAlignmentsInReadGroup
    File? outputAlignmentsInFileReadGroups
    Int? mapqThreshold
    String? outputAlignmentsInLibrary
    Int? outputAlignmentsMeetingCIGARThreshold
    String? outputAlignmentsWithBitsSet
    String? doNotOutputAlignmentsWithBitsSet
    String? doNotOutputAlignmentsWithAllBitsSet
    String? readTagToExclude
    Boolean? collapseBackwardCIGAR
    Float? subsamplingProportion
    Int? threads
    File sam
    File? reference
    File? reference_fai
    File? reference_amb
    File? reference_ann
    File? reference_bwt
    File? reference_pac
    File? reference_sa
    File? reference_dict
    String? outputFilename
    Array[String]? regions
  }
  command <<<
    set -e
    samtools view \
      '-S' \
      '-h' \
      '-b' \
      ~{if (defined(cramOutput) && select_first([cramOutput])) then "-C" else ""} \
      ~{if (defined(compressedBam) && select_first([compressedBam])) then "-1" else ""} \
      ~{if (defined(uncompressedBam) && select_first([uncompressedBam])) then "-u" else ""} \
      ~{if (defined(onlyOutputHeader) && select_first([onlyOutputHeader])) then "-H" else ""} \
      ~{if (defined(countAlignments) && select_first([countAlignments])) then "-c" else ""} \
      ~{if defined(writeAlignments) then ("-U '" + writeAlignments + "'") else ""} \
      ~{if defined(inputTSV) then ("-t '" + inputTSV + "'") else ""} \
      ~{if defined(onlyOverlapping) then ("-L '" + onlyOverlapping + "'") else ""} \
      ~{if (defined(useMultiRegionIterator) && select_first([useMultiRegionIterator])) then "-M" else ""} \
      ~{if defined(outputAlignmentsInReadGroup) then ("-r '" + outputAlignmentsInReadGroup + "'") else ""} \
      ~{if defined(outputAlignmentsInFileReadGroups) then ("-R '" + outputAlignmentsInFileReadGroups + "'") else ""} \
      ~{if defined(mapqThreshold) then ("-q " + mapqThreshold) else ''} \
      ~{if defined(outputAlignmentsInLibrary) then ("-l '" + outputAlignmentsInLibrary + "'") else ""} \
      ~{if defined(outputAlignmentsMeetingCIGARThreshold) then ("-m " + outputAlignmentsMeetingCIGARThreshold) else ''} \
      ~{if defined(outputAlignmentsWithBitsSet) then ("-f '" + outputAlignmentsWithBitsSet + "'") else ""} \
      ~{if defined(doNotOutputAlignmentsWithBitsSet) then ("-F '" + doNotOutputAlignmentsWithBitsSet + "'") else ""} \
      ~{if defined(doNotOutputAlignmentsWithAllBitsSet) then ("-G '" + doNotOutputAlignmentsWithAllBitsSet + "'") else ""} \
      ~{if defined(readTagToExclude) then ("-x '" + readTagToExclude + "'") else ""} \
      ~{if (defined(collapseBackwardCIGAR) && select_first([collapseBackwardCIGAR])) then "-B" else ""} \
      ~{if defined(subsamplingProportion) then ("-s " + subsamplingProportion) else ''} \
      ~{if defined(threads) then ("-@ " + threads) else ''} \
      -o '~{select_first([outputFilename, "generated.bam"])}' \
      ~{if defined(reference) then ("-T '" + reference + "'") else ""} \
      ~{sam} \
      ~{if (defined(regions) && length(select_first([regions])) > 0) then "'" + sep("' '", select_first([regions])) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/samtools:1.9--h8571acd_11"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.bam"])
  }
}