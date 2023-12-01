
# trim read 1 and read 2 adapter sequeunce with cutadapt
task TrimAdapters {
  input {
    File fastq_input_read1
    File fastq_input_read2
    Int min_length
    Int quality_cutoff
    String adapter_seq_read1
    String adapter_seq_read2
    String output_base_name
    String docker_image = "quay.io/broadinstitute/cutadapt:1.18"
  }

  parameter_meta {
    fastq_input_read1: "read 1 fastq file as input for the pipeline"
    fastq_input_read2: "read 2 fastq file as input for the pipeline"
    min_length: "the minimum legnth for trimming. Reads that are too short even before adapter removal are also discarded"
    quality_cutoff: "cutadapt option to trim low-quality ends from reads before adapter removal"
    adapter_seq_read1: "cutadapt option for the sequence adapter for read 1 fastq"
    adapter_seq_read2: "cutadapt option for the sequence adapter for read 2 fastq"
    output_base_name: "base name to be used for the output of the task"
    docker_image: "the docker image using cutadapt to be used (default: quay.io/broadinstitute/cutadapt:1.18)"
  }

  # runtime requirements based upon input file size
  Float input_size = size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB")
  Int disk_size = ceil(2 * (if input_size < 1 then 1 else input_size))

  # output names for trimmed reads
  String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
  String fastq_trimmed_adapter_output_name_read2 = output_base_name + ".R2.trimmed_adapters.fastq.gz"

  # using cutadapt to trim off sequence adapters
  command {
    set -euo pipefail

    # fastq's, "-f", -A for paired adapters read 2"
    cutadapt \
      -f fastq \
      --minimum-length ~{min_length} \
      --quality-cutoff ~{quality_cutoff} \
      --adapter ~{adapter_seq_read1} \
      -A ~{adapter_seq_read2} \
      --output ~{fastq_trimmed_adapter_output_name_read1} \
      --paired-output ~{fastq_trimmed_adapter_output_name_read2} \
      ~{fastq_input_read1} ~{fastq_input_read2}
  }

  # use docker image for given tool cutadapat
  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File fastq_trimmed_adapter_output_read1 = fastq_trimmed_adapter_output_name_read1
    File fastq_trimmed_adapter_output_read2 = fastq_trimmed_adapter_output_name_read2
    File? monitoring_log = "monitoring.log"
  }
}
