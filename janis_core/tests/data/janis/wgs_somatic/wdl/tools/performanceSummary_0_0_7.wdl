version development

task performanceSummary {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File flagstat
    File collectInsertSizeMetrics
    File coverage
    String? outputPrefix
    File? targetFlagstat
    File? rmdupFlagstat
    Boolean? genome
  }
  command <<<
    set -e
    performance_summary.py \
      --flagstat '~{flagstat}' \
      --collect_insert_metrics '~{collectInsertSizeMetrics}' \
      --coverage '~{coverage}' \
      -o '~{select_first([outputPrefix, "generated.csv"])}' \
      ~{if defined(targetFlagstat) then ("--target_flagstat '" + targetFlagstat + "'") else ""} \
      ~{if defined(rmdupFlagstat) then ("--rmdup_flagstat '" + rmdupFlagstat + "'") else ""} \
      ~{if (defined(genome) && select_first([genome])) then "--genome" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/pmacutil:0.0.7"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = (select_first([outputPrefix, "generated.csv"]) + ".csv")
  }
}