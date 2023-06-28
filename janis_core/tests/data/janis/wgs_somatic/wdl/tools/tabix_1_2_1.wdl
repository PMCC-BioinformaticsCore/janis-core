version development

task tabix {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File inp
    String? preset
    Boolean? zeroBased
    Int? begin
    String? comment
    Boolean? csi
    Int? end
    Boolean? force
    Int? minShift
    Int? sequence
    Int? skipLines
    Boolean? printHeader
    Boolean? onlyHeader
    Boolean? listChroms
    File? reheader
    File? regions
    File? targets
  }
  command <<<
    set -e
    cp -f '~{inp}' '.'
    tabix \
      ~{if (defined(zeroBased) && select_first([zeroBased])) then "--zero-based" else ""} \
      ~{if (defined(csi) && select_first([csi])) then "--csi" else ""} \
      ~{if (defined(force) && select_first([force])) then "--force" else ""} \
      ~{if defined(minShift) then ("--min-shift " + minShift) else ''} \
      ~{if (defined(printHeader) && select_first([printHeader])) then "--print-header" else ""} \
      ~{if (defined(onlyHeader) && select_first([onlyHeader])) then "--only-header" else ""} \
      ~{if (defined(listChroms) && select_first([listChroms])) then "--list-chroms" else ""} \
      ~{if defined(reheader) then ("--reheader '" + reheader + "'") else ""} \
      ~{if defined(select_first([preset, "vcf"])) then ("--preset '" + select_first([preset, "vcf"]) + "'") else ""} \
      ~{if defined(sequence) then ("--sequence " + sequence) else ''} \
      ~{if defined(begin) then ("--begin " + begin) else ''} \
      ~{if defined(end) then ("--end " + end) else ''} \
      ~{if defined(skipLines) then ("--skip-lines " + skipLines) else ''} \
      ~{if defined(comment) then ("--comment '" + comment + "'") else ""} \
      '~{basename(inp)}' \
      ~{if defined(regions) then ("--regions '" + regions + "'") else ""} \
      ~{if defined(targets) then ("--targets '" + targets + "'") else ""}
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
    File out = basename(inp)
    File out_tbi = basename(inp) + ".tbi"
  }
}