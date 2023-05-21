version development

task bcftoolsAnnotate {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File vcf
    String? outputFilename
    File? annotations
    String? collapse
    Array[String]? columns
    String? exclude
    File? headerLines
    String? setId
    String? include
    Boolean? keepSites
    String? markSites
    String? outputType
    String? regions
    File? regionsFile
    File? renameChrs
    Array[File]? samples
    File? samplesFile
    Int? threads
    Array[String]? remove
  }
  command <<<
    set -e
    bcftools annotate \
      --output '~{select_first([outputFilename, "generated.vcf"])}' \
      ~{if defined(annotations) then ("--annotations '" + annotations + "'") else ""} \
      ~{if defined(collapse) then ("--collapse '" + collapse + "'") else ""} \
      ~{if (defined(columns) && length(select_first([columns])) > 0) then "--columns '" + sep("' '", select_first([columns])) + "'" else ""} \
      ~{if defined(exclude) then ("--exclude '" + exclude + "'") else ""} \
      ~{if defined(headerLines) then ("--header-lines '" + headerLines + "'") else ""} \
      ~{if defined(setId) then ("--set-id '" + setId + "'") else ""} \
      ~{if defined(include) then ("--include '" + include + "'") else ""} \
      ~{if (defined(keepSites) && select_first([keepSites])) then "--keep-sites" else ""} \
      ~{if defined(markSites) then ("--mark-sites '" + markSites + "'") else ""} \
      ~{if defined(outputType) then ("--output-type '" + outputType + "'") else ""} \
      ~{if defined(regions) then ("--regions '" + regions + "'") else ""} \
      ~{if defined(regionsFile) then ("--regions-file '" + regionsFile + "'") else ""} \
      ~{if defined(renameChrs) then ("--rename-chrs '" + renameChrs + "'") else ""} \
      ~{if (defined(samples) && length(select_first([samples])) > 0) then "--samples '" + sep("' '", select_first([samples])) + "'" else ""} \
      ~{if defined(samplesFile) then ("--samples-file '" + samplesFile + "'") else ""} \
      ~{if defined(threads) then ("--threads " + threads) else ''} \
      ~{if (defined(remove) && length(select_first([remove])) > 0) then "--remove '" + sep("' '", select_first([remove])) + "'" else ""} \
      '~{vcf}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biocontainers/bcftools:v1.5_cv2"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.vcf"])
  }
}