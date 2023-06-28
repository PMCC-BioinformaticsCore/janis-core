version development

task manta {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File? config
    File bam
    File bam_bai
    String? runDir
    File reference
    File reference_fai
    File? tumorBam
    File? tumorBam_bai
    Boolean? exome
    Boolean? rna
    Boolean? unstrandedRNA
    Boolean? outputContig
    File? callRegions
    File? callRegions_tbi
    String? mode
    Boolean? quiet
    String? queue
    Int? memgb
    String? maxTaskRuntime
  }
  command <<<
    set -e
     \
      configManta.py \
      ~{if defined(config) then ("--config " + config) else ''} \
      --bam ~{bam} \
      --runDir ~{select_first([runDir, "generated"])} \
      --referenceFasta ~{reference} \
      ~{if defined(tumorBam) then ("--tumorBam " + tumorBam) else ''} \
      ~{if (defined(exome) && select_first([exome])) then "--exome" else ""} \
      ~{if (defined(rna) && select_first([rna])) then "--rna" else ""} \
      ~{if (defined(unstrandedRNA) && select_first([unstrandedRNA])) then "--unstrandedRNA" else ""} \
      ~{if (defined(outputContig) && select_first([outputContig])) then "--outputContig" else ""} \
      ~{if defined(callRegions) then ("--callRegions " + callRegions) else ''} \
      ;~{select_first([runDir, "generated"])}/runWorkflow.py \
      ~{if defined(select_first([mode, "local"])) then ("--mode " + select_first([mode, "local"])) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      ~{if defined(queue) then ("--queue " + queue) else ''} \
      ~{if defined(memgb) then ("--memGb " + memgb) else ''} \
      ~{if defined(maxTaskRuntime) then ("--maxTaskRuntime " + maxTaskRuntime) else ''} \
      -j ~{select_first([runtime_cpu, 4, 1])}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 4, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/manta:1.5.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4, 4])}G"
    preemptible: 2
  }
  output {
    File python = (select_first([runDir, "generated"]) + "/runWorkflow.py")
    File pickle = (select_first([runDir, "generated"]) + "/runWorkflow.py.config.pickle")
    File candidateSV = (select_first([runDir, "generated"]) + "/results/variants/candidateSV.vcf.gz")
    File candidateSV_tbi = (select_first([runDir, "generated"]) + "/results/variants/candidateSV.vcf.gz") + ".tbi"
    File candidateSmallIndels = (select_first([runDir, "generated"]) + "/results/variants/candidateSmallIndels.vcf.gz")
    File candidateSmallIndels_tbi = (select_first([runDir, "generated"]) + "/results/variants/candidateSmallIndels.vcf.gz") + ".tbi"
    File diploidSV = (select_first([runDir, "generated"]) + "/results/variants/diploidSV.vcf.gz")
    File diploidSV_tbi = (select_first([runDir, "generated"]) + "/results/variants/diploidSV.vcf.gz") + ".tbi"
    File alignmentStatsSummary = (select_first([runDir, "generated"]) + "/results/stats/alignmentStatsSummary.txt")
    File svCandidateGenerationStats = (select_first([runDir, "generated"]) + "/results/stats/svCandidateGenerationStats.tsv")
    File svLocusGraphStats = (select_first([runDir, "generated"]) + "/results/stats/svLocusGraphStats.tsv")
    File? somaticSVs = (select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz")
    File? somaticSVs_tbi = if defined((select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz")) then ((select_first([runDir, "generated"]) + "/results/variants/somaticSV.vcf.gz") + ".tbi") else None
  }
}