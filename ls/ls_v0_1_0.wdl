version development

task ls {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disk
    File inp
    File inp_sec
  }

  command <<<
    set -e
    cp -f '~{inp_sec}' $(echo '~{inp}' | sed 's/\.[^.]*$//').sec
    ls
    if [ -f $(echo '~{inp}' | sed 's/\.[^.]*$//').sec ]; then ln -f $(echo '~{inp}' | sed 's/\.[^.]*$//').sec $(echo '~{inp}' ).sec; fi
  >>>

  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
    docker: "ubuntu:latest"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }

  output {
    File std = stdout()
    File out = inp
    File out_sec = inp + ".sec"
  }

}