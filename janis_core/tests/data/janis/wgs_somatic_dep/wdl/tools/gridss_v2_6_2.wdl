version development

task gridss {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] bams
    Array[File] bams_bai
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputFilename
    String? assemblyFilename
    Int? threads
    File? blacklist
    String? tmpdir
  }
  command <<<
    set -e
    /opt/gridss/gridss.sh \
      ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
      ~{if defined(select_first([tmpdir, "./TMP"])) then ("--workingdir '" + select_first([tmpdir, "./TMP"]) + "'") else ""} \
      --reference '~{reference}' \
      --output '~{select_first([outputFilename, "generated.svs.vcf"])}' \
      --assembly '~{select_first([assemblyFilename, "generated.assembled.bam"])}' \
      ~{if defined(blacklist) then ("--blacklist '" + blacklist + "'") else ""} \
      ~{if length(bams) > 0 then "'" + sep("' '", bams) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 8, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "gridss/gridss:2.6.2"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 31, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.svs.vcf"])
    File assembly = select_first([assemblyFilename, "generated.assembled.bam"])
  }
}