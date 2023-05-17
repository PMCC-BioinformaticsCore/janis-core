version development

task VcfTools {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File vcf
    String? outputFilename
    Boolean? removeFileteredAll
    Boolean? recode
    Boolean? recodeINFOAll
  }
  command <<<
    set -e
     vcftools \
      --vcf '~{vcf}' \
      --out '~{select_first([outputFilename, "generated"])}' \
      ~{if (defined(removeFileteredAll) && select_first([removeFileteredAll])) then "--remove-filtered-all" else ""} \
      ~{if (defined(recode) && select_first([recode])) then "--recode" else ""} \
      ~{if (defined(recodeINFOAll) && select_first([recodeINFOAll])) then "--recode-INFO-all" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biocontainers/vcftools:v0.1.16-1-deb_cv1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = (select_first([outputFilename, "generated"]) + ".recode.vcf")
  }
}