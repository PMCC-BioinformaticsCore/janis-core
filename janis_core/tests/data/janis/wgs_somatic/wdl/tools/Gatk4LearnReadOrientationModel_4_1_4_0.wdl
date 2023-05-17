version development

task Gatk4LearnReadOrientationModel {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    Array[File] f1r2CountsFiles
    Int? numEmIterations
    String? modelFileOut
  }
  command <<<
    set -e
    gatk LearnReadOrientationModel \
      --java-options '-Xmx~{((select_first([runtime_memory, 32, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if length(f1r2CountsFiles) > 0 then "-I '" + sep("' -I '", f1r2CountsFiles) + "'" else ""} \
      ~{if defined(select_first([numEmIterations, 30])) then ("--num-em-iterations " + select_first([numEmIterations, 30])) else ''} \
      -O '~{select_first([modelFileOut, "generated.tar.gz"])}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.4.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 32, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([modelFileOut, "generated.tar.gz"])
  }
}