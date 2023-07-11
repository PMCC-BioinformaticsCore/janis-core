version development

task FilterVardictSomaticVcf {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File? vcf
    String? outputFilename
  }
  command <<<
    set -e
     \
      bcftools filter -e 'STATUS=\"GERMLINE\"' -o - \
      ~{if defined(vcf) then ("'" + vcf + "'") else ""} \
      | bcftools filter -i 'FILTER==\"PASS\"' \
      -o ~{select_first([outputFilename, "generated.filter.vcf"])}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.filter.vcf"])
  }
}