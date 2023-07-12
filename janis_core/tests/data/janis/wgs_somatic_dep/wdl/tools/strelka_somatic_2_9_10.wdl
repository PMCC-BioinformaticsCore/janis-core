version development

task strelka_somatic {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File normalBam
    File normalBam_bai
    File tumorBam
    File tumorBam_bai
    File reference
    File reference_fai
    String? rundir
    Array[String]? region
    File? config
    Boolean? outputcallableregions
    Array[File]? indelCandidates
    Array[File]? indelCandidates_tbi
    Array[File]? forcedgt
    Array[File]? forcedgt_tbi
    Boolean? targeted
    Boolean? exome
    File? callRegions
    File? callRegions_tbi
    File? noisevcf
    File? noisevcf_tbi
    Int? scansizemb
    Int? callmemmb
    Boolean? retaintempfiles
    Boolean? disableevs
    Boolean? reportevsfeatures
    File? snvscoringmodelfile
    File? indelscoringmodelfile
    String? mode
    String? queue
    String? memGb
    Boolean? quiet
  }
  command <<<
    set -e
     \
      'configureStrelkaSomaticWorkflow.py' \
      --normalBam='~{normalBam}' \
      --tumourBam='~{tumorBam}' \
      --referenceFasta='~{reference}' \
      --runDir='~{select_first([rundir, "generated"])}' \
      ~{if (defined(region) && length(select_first([region])) > 0) then "--region '" + sep("' --region '", select_first([region])) + "'" else ""} \
      ~{if defined(config) then ("--config='" + config + "'") else ""} \
      ~{if (defined(outputcallableregions) && select_first([outputcallableregions])) then "--outputCallableRegions" else ""} \
      ~{if (defined(indelCandidates) && length(select_first([indelCandidates])) > 0) then "--indelCandidates='" + sep("' --indelCandidates='", select_first([indelCandidates])) + "'" else ""} \
      ~{if (defined(forcedgt) && length(select_first([forcedgt])) > 0) then "--forcedGT='" + sep("' --forcedGT='", select_first([forcedgt])) + "'" else ""} \
      ~{if (defined(targeted) && select_first([targeted])) then "--targeted" else ""} \
      ~{if (defined(exome) && select_first([exome])) then "--exome" else ""} \
      ~{if defined(callRegions) then ("--callRegions='" + callRegions + "'") else ""} \
      ~{if defined(noisevcf) then ("--noiseVcf='" + noisevcf + "'") else ""} \
      ~{if defined(scansizemb) then ("--scanSizeMb=" + scansizemb) else ''} \
      ~{if defined(callmemmb) then ("--callMemMb=" + callmemmb) else ''} \
      ~{if select_first([retaintempfiles, false]) then "--retainTempFiles" else ""} \
      ~{if (defined(disableevs) && select_first([disableevs])) then "--disableEVS" else ""} \
      ~{if (defined(reportevsfeatures) && select_first([reportevsfeatures])) then "--reportEVSFeatures" else ""} \
      ~{if defined(snvscoringmodelfile) then ("--snvScoringModelFile='" + snvscoringmodelfile + "'") else ""} \
      ~{if defined(indelscoringmodelfile) then ("--indelScoringModelFile='" + indelscoringmodelfile + "'") else ""} \
      ;~{select_first([rundir, "generated"])}/runWorkflow.py \
      ~{if defined(select_first([mode, "local"])) then ("--mode " + select_first([mode, "local"])) else ''} \
      ~{if defined(queue) then ("--queue " + queue) else ''} \
      ~{if defined(memGb) then ("--memGb " + memGb) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      --jobs ~{select_first([runtime_cpu, 4, 1])}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 4, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/strelka:2.9.10"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4, 4])}G"
    preemptible: 2
  }
  output {
    File configPickle = (select_first([rundir, "generated"]) + "/runWorkflow.py.config.pickle")
    File script = (select_first([rundir, "generated"]) + "/runWorkflow.py")
    File stats = (select_first([rundir, "generated"]) + "/results/stats/runStats.tsv")
    File indels = (select_first([rundir, "generated"]) + "/results/variants/somatic.indels.vcf.gz")
    File indels_tbi = (select_first([rundir, "generated"]) + "/results/variants/somatic.indels.vcf.gz") + ".tbi"
    File snvs = (select_first([rundir, "generated"]) + "/results/variants/somatic.snvs.vcf.gz")
    File snvs_tbi = (select_first([rundir, "generated"]) + "/results/variants/somatic.snvs.vcf.gz") + ".tbi"
  }
}