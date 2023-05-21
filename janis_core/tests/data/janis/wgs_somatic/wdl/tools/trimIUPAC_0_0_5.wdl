version development

task trimIUPAC {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File vcf
    String? outputFilename
  }
  command <<<
    set -e
    trimIUPAC.py \
      '~{vcf}' \
      '~{select_first([outputFilename, "generated.trimmed.vcf"])}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/pmacutil:0.0.5"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 1, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.trimmed.vcf"])
  }
}