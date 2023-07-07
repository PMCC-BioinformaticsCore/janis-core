version development

task extractStrelkaSomaticADDP {
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
    extract_strelka_somatic_DP_AF.py \
      -i '~{vcf}' \
      -o '~{select_first([outputFilename, "generated.vcf"])}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/pmacutil:0.1.1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.vcf"])
  }
}