version development

task bcftoolssort {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File vcf
    String? outputFilename
    String? outputType
    String? tempDir
  }
  command <<<
    set -e
    bcftools sort \
      --output-file '~{select_first([outputFilename, "generated.sorted.vcf.gz"])}' \
      ~{if defined(select_first([outputType, "z"])) then ("--output-type '" + select_first([outputType, "z"]) + "'") else ""} \
      ~{if defined(tempDir) then ("--temp-dir '" + tempDir + "'") else ""} \
      ~{vcf}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.sorted.vcf.gz"])
  }
}