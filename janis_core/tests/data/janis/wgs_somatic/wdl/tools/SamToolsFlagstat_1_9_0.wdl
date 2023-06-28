version development

task SamToolsFlagstat {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File bam
    Int? threads
  }
  command <<<
    set -e
    samtools flagstat \
      ~{if defined(threads) then ("-@ " + threads) else ''} \
      '~{bam}'
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
    File out = stdout()
  }
}