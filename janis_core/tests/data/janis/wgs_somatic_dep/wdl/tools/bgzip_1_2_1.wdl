version development

task bgzip {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File file
    String? outputFilename
    Int? offset
    Boolean? stdout
    Boolean? decompress
    Boolean? force
    Boolean? help
    Boolean? index
    File? indexName
    Int? compress
    Boolean? reindex
    Boolean? rebgzip
    Int? size
    Int? threads
  }
  command <<<
    set -e
    bgzip \
      ~{if defined(offset) then ("--offset " + offset) else ''} \
      ~{if select_first([stdout, true]) then "--stdout" else ""} \
      ~{if (defined(decompress) && select_first([decompress])) then "--decompress" else ""} \
      ~{if (defined(force) && select_first([force])) then "--force" else ""} \
      ~{if (defined(help) && select_first([help])) then "--help" else ""} \
      ~{if (defined(index) && select_first([index])) then "--index" else ""} \
      ~{if defined(indexName) then ("--index-name '" + indexName + "'") else ""} \
      ~{if defined(compress) then ("--compress " + compress) else ''} \
      ~{if (defined(reindex) && select_first([reindex])) then "--reindex" else ""} \
      ~{if (defined(rebgzip) && select_first([rebgzip])) then "--rebgzip" else ""} \
      ~{if defined(size) then ("--size " + size) else ''} \
      ~{if defined(threads) then ("--threads " + threads) else ''} \
      '~{file}' \
      > \
      ~{select_first([outputFilename, "~{basename(file)}.gz"])}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "biodckrdev/htslib:1.2.1"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "~{basename(file)}.gz"])
  }
}