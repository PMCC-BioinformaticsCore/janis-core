version development

task UncompressArchive {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File file
    Boolean? stdout
    Boolean? decompress
    Boolean? force
    Boolean? keep
    Boolean? list
    Boolean? noName
    Boolean? name
    Boolean? quiet
    Boolean? recursive
    String? suffix
    Boolean? test
    Boolean? fast
    Boolean? best
    Boolean? rsyncable
  }
  command <<<
    set -e
    gunzip \
      ~{if select_first([stdout, true]) then "-c" else ""} \
      ~{if (defined(decompress) && select_first([decompress])) then "-d" else ""} \
      ~{if (defined(force) && select_first([force])) then "-f" else ""} \
      ~{if (defined(keep) && select_first([keep])) then "-k" else ""} \
      ~{if (defined(list) && select_first([list])) then "-l" else ""} \
      ~{if (defined(noName) && select_first([noName])) then "-n" else ""} \
      ~{if (defined(name) && select_first([name])) then "-N" else ""} \
      ~{if (defined(quiet) && select_first([quiet])) then "-q" else ""} \
      ~{if (defined(recursive) && select_first([recursive])) then "-r" else ""} \
      ~{if defined(suffix) then ("-s '" + suffix + "'") else ""} \
      ~{if (defined(test) && select_first([test])) then "-t" else ""} \
      ~{if (defined(fast) && select_first([fast])) then "-1" else ""} \
      ~{if (defined(best) && select_first([best])) then "-9" else ""} \
      ~{if (defined(rsyncable) && select_first([rsyncable])) then "--rsyncable" else ""} \
      '~{file}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "ubuntu:latest"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}