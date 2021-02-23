version development

task cat {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File? file
    Array[File]? files
    Boolean? number_output
    Boolean? number_non_blank
    Boolean? disable_output_buffer
    Boolean? squeeze
    Boolean? display_nonprint_and_eol_chars
    Boolean? display_nonprint_and_tab_chars
    Boolean? display_nonprint_chars
  }
  command <<<
    set -e
    cat \
      ~{if (defined(number_output) && select_first([number_output])) then "-n" else ""} \
      ~{if (defined(number_non_blank) && select_first([number_non_blank])) then "-b" else ""} \
      ~{if (defined(disable_output_buffer) && select_first([disable_output_buffer])) then "-u" else ""} \
      ~{if (defined(squeeze) && select_first([squeeze])) then "-s" else ""} \
      ~{if (defined(display_nonprint_and_eol_chars) && select_first([display_nonprint_and_eol_chars])) then "-e" else ""} \
      ~{if (defined(display_nonprint_and_tab_chars) && select_first([display_nonprint_and_tab_chars])) then "-t" else ""} \
      ~{if (defined(display_nonprint_chars) && select_first([display_nonprint_chars])) then "-v" else ""} \
      ~{if (defined(files) && length(select_first([files])) > 0) then "'" + sep("' '", select_first([files])) + "'" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}