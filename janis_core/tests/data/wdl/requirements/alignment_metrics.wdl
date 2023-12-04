version 1.0

task alignment_metrics {
  input {
    File   aligned_bam
    File   ref_fasta
    File?  primers_bed

    Int?   machine_mem_gb
    String docker = "quay.io/broadinstitute/viral-core:2.1.33"
  }

  String out_basename = basename(aligned_bam, ".bam")

  command <<<
    set -e
    MEM_MB=$(free -m | head -2 | tail -1 | awk '{print $4}')
    XMX=$(echo "-Xmx"$MEM_MB"m")
    echo "Requesting $MEM_MB MB of RAM for Java"

    # requisite Picard fasta indexing
    cp "~{ref_fasta}" reference.fasta
    picard $XMX CreateSequenceDictionary -R reference.fasta

    # get Picard metrics and clean up the junky outputs
    picard $XMX CollectRawWgsMetrics \
      -R reference.fasta \
      -I "~{aligned_bam}" \
      -O picard_raw.raw_wgs_metrics.txt
    grep -v \# picard_raw.raw_wgs_metrics.txt | grep . | head -2 > picard_clean.raw_wgs_metrics.txt

    picard $XMX CollectAlignmentSummaryMetrics \
      -R reference.fasta \
      -I "~{aligned_bam}" \
      -O picard_raw.alignment_metrics.txt
    grep -v \# picard_raw.alignment_metrics.txt | grep . | head -4 > picard_clean.alignment_metrics.txt 

    picard $XMX CollectInsertSizeMetrics \
      -I "~{aligned_bam}" \
      -O picard_raw.insert_size_metrics.txt \
      -H picard_raw.insert_size_metrics.pdf \
      --INCLUDE_DUPLICATES true
    grep -v \# picard_raw.insert_size_metrics.txt | grep . | head -2 > picard_clean.insert_size_metrics.txt

    # prepend the sample name in order to facilitate tsv joining later
    SAMPLE=$(samtools view -H "~{aligned_bam}" | grep ^@RG | perl -lape 's/^@RG.*SM:(\S+).*$/$1/' | sort | uniq)
    echo -e "sample_sanitized\tbam" > prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.raw_wgs_metrics.txt > "~{out_basename}".raw_wgs_metrics.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.alignment_metrics.txt > "~{out_basename}".alignment_metrics.txt
    echo -e "sample_sanitized\tbam" > prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.insert_size_metrics.txt > "~{out_basename}".insert_size_metrics.txt

    # actually don't know how to do CollectTargetedPcrMetrics yet
    if [ -n "~{primers_bed}" ]; then
      picard $XMX BedToIntervalList \
        -I "~{primers_bed}" \
        -O primers.interval.list \
        -SD reference.dict
    fi
  >>>

  output {
    File wgs_metrics         = "~{out_basename}.raw_wgs_metrics.txt"
    File aln_metrics         = "~{out_basename}.alignment_metrics.txt"
    File insert_size_metrics = "~{out_basename}.insert_size_metrics.txt"
  }

  runtime {
    docker: "~{docker}"
    memory: select_first([machine_mem_gb, 13]) + " GB"
    cpu: 2
    disks: "local-disk 150 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}