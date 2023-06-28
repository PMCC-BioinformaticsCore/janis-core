version development

task SplitMultiAllele {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File vcf
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputFilename
  }
  command <<<
    set -e
     \
      vt decompose -s \
      ~{vcf} \
      | vt normalize -n -q - \
      -r ~{reference} \
      -o ~{select_first([outputFilename, "generated.norm.vcf"])}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "heuermh/vt"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.norm.vcf"])
  }
}