version development

task addBamStats {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File? mpileup
    File? normalMpileup
    File? tumorMpileup
    String? normalID
    String? tumorID
    File inputVcf
    String? outputFilename
    String type
  }
  command <<<
    set -e
    add_bam_stats.py \
      ~{if defined(mpileup) then ("--mpileup '" + mpileup + "'") else ""} \
      ~{if defined(normalMpileup) then ("--normal_mpileup '" + normalMpileup + "'") else ""} \
      ~{if defined(tumorMpileup) then ("--tumor_mpileup '" + tumorMpileup + "'") else ""} \
      ~{if defined(normalID) then ("--normal_id '" + normalID + "'") else ""} \
      ~{if defined(tumorID) then ("--tumor_id '" + tumorID + "'") else ""} \
      -i '~{inputVcf}' \
      -o '~{select_first([outputFilename, "generated.addbamstats.vcf"])}' \
      --type '~{type}'
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
    File out = select_first([outputFilename, "generated.addbamstats.vcf"])
  }
}