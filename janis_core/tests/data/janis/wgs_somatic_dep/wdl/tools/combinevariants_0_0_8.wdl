version development

task combinevariants {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    String? outputFilename
    Array[File] vcfs
    String type
    Array[String]? columns
    String? normal
    String? tumor
    Int? priority
  }
  command <<<
    set -e
    combine_vcf.py \
      -o '~{select_first([outputFilename, "generated.combined.vcf"])}' \
      ~{if length(vcfs) > 0 then "-i '" + sep("' -i '", vcfs) + "'" else ""} \
      --type '~{type}' \
      ~{if (defined(columns) && length(select_first([columns])) > 0) then "--columns '" + sep("','", select_first([columns])) + "'" else ""} \
      ~{if defined(normal) then ("--normal '" + normal + "'") else ""} \
      ~{if defined(tumor) then ("--tumor '" + tumor + "'") else ""} \
      ~{if defined(priority) then ("--priority " + priority) else ''}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/pmacutil:0.0.8"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.combined.vcf"])
  }
}