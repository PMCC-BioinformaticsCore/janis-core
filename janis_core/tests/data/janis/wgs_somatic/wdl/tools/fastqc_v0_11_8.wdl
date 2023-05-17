version development

task fastqc {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] reads
    String? outdir
    Boolean? casava
    Boolean? nano
    Boolean? nofilter
    Boolean? extract
    String? java
    Boolean? noextract
    Boolean? nogroup
    String? format
    Int? threads
    File? contaminants
    File? adapters
    File? limits
    Int? kmers
    Boolean? quiet
    String? dir
  }
  command <<<
    set -e
    fastqc \
      ~{if defined(select_first([outdir, "."])) then ("--outdir '" + select_first([outdir, "."]) + "'") else ""} \
      ~{if (defined(casava) && select_first([casava])) then "--casava" else ""} \
      ~{if (defined(nano) && select_first([nano])) then "--nano" else ""} \
      ~{if (defined(nofilter) && select_first([nofilter])) then "--nofilter" else ""} \
      ~{if select_first([extract, true]) then "--extract" else ""} \
      ~{if defined(java) then ("--java '" + java + "'") else ""} \
      ~{if (defined(noextract) && select_first([noextract])) then "--noextract" else ""} \
      ~{if (defined(nogroup) && select_first([nogroup])) then "--nogroup" else ""} \
      ~{if defined(format) then ("--format '" + format + "'") else ""} \
      ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
      ~{if defined(contaminants) then ("--contaminants '" + contaminants + "'") else ""} \
      ~{if defined(adapters) then ("--adapters '" + adapters + "'") else ""} \
      ~{if defined(limits) then ("--limits '" + limits + "'") else ""} \
      ~{if defined(kmers) then ("--kmers " + kmers) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      ~{if defined(dir) then ("--dir '" + dir + "'") else ""} \
      ~{if length(reads) > 0 then "'" + sep("' '", reads) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/fastqc:0.11.8--2"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    Array[File] out = glob("*.zip")
    Array[File] datafile = glob("*/fastqc_data.txt")
  }
}