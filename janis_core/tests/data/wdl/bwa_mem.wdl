version 1.0

# align the two trimmed fastq as piared end data using BWA
task BWAPairedEndAlignment {
  input {
    File fastq_input_read1
    File fastq_input_read2
    File tar_bwa_reference
    String read_group_id
    String read_group_sample_name
    Int cpu
    String output_base_name
    String docker_image = "quay.io/humancellatlas/snaptools:0.0.1"
  }

  parameter_meta {
    fastq_input_read1: "the trimmed read 1 fastq file as input for the aligner"
    fastq_input_read2: "the trimmed read 1 fastq file as input for the aligner"
    tar_bwa_reference: "the pre built tar file containing the reference fasta and cooresponding reference files for the BWA aligner"
    read_group_id: "the read group id to be added upon alignment"
    read_group_sample_name: "the read group sample to be added upon alignment"
    cpu: "the number of cpu cores to use during alignment"
    output_base_name: "basename to be used for the output of the task"
    docker_image: "the docker image using BWA to be used (default: quay.io/humancellatlas/snaptools:0.0.1)"
  }

  # runtime requirements based upon input file size
  Float input_size = size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB") + size(tar_bwa_reference, "GiB")
  Int disk_size = ceil(3.25 * (if input_size < 1 then 1 else input_size))

  String sam_aligned_output_name = output_base_name + ".aligned.sam"

  # sort with samtools
  command {
    set -euo pipefail

    # prepare reference
    declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
    tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
    rm "~{tar_bwa_reference}" || /bin/true

    # align w/ BWA: -t for number of cores
    bwa \
      mem \
      -R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" \ 
      -t ~{cpu} \
      $REF_DIR/genome.fa \
      <(zcat ~{fastq_input_read1}) <(zcat ~{fastq_input_read2}) \
      > ~{sam_aligned_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    memory: "3.75 GiB"
  }

  output {
    File sam_aligned_output = sam_aligned_output_name
    File? monitoring_log = "monitoring.log"
  }
}