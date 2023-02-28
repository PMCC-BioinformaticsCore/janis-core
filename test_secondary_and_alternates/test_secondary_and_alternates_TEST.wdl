version development

task test_secondary_and_alternates {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disk
    File inp
    File inp_file
  }

  command <<<
    set -e
    cp -f '~{inp}' '.'
    cp -f '~{inp_file}' .
    cat \
      '~{basename(inp)}'
  >>>

  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
    docker: "ubtunu"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }

  output {
    File out = basename(inp)
    File out_file = sub(sub(basename(inp), "\\.txt$", ".file"), "\\.text$", ".file")
  }

}