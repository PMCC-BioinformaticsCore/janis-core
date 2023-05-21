nextflow.enable.dsl=2

process MANTA {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/manta:1.4.0"
    publishDir "${params.outdir}/run_manta/manta"
    cpus "${params.run_manta.manta.cpus}"
    memory "${params.run_manta.manta.memory}"

    input:
    tuple path(primary), path(tbi)
    tuple path(primary), path(crai)
    tuple path(primary), path(crai)
    tuple path(primary), path(dict), path(fai)
    val output_basename
    val cores
    val ram

    output:
    tuple path("*SV.vcf.gz"), path("*.tbi"), emit: output_sv
    tuple path("*SmallIndels.vcf.gz"), path("*.tbi"), emit: small_indels

    script:
    """
    /manta-1.4.0.centos6_x86_64/bin/configManta.py \
    <js>${  var std = " --ref " + inputs.reference.path + " --callRegions " + inputs.hg38_strelka_bed.path + " --runDir=./ && ./runWorkflow.py -m local -j " + inputs.cores + " ";  var mv = " && mv results/variants/";  if (typeof inputs.input_tumor_cram === 'undefined' || inputs.input_tumor_cram === null){    var mv_cmd = mv + "diploidSV.vcf.gz " +  inputs.output_basename + ".manta.diploidSV.vcf.gz" + mv + "diploidSV.vcf.gz.tbi " + inputs.output_basename + ".manta.diploidSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";    return "--bam ".concat(inputs.input_normal_cram.path, std, mv_cmd);  }  else if (typeof inputs.input_normal_cram === 'undefined' || inputs.input_normal_cram === null){    var mv_cmd = mv + "tumorSV.vcf.gz " + inputs.output_basename + ".manta.tumorSV.vcf.gz" + mv + "tumorSV.vcf.gz.tbi " + inputs.output_basename + ".manta.tumorSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";    return "--tumorBam " + inputs.input_tumor_cram.path + std + mv_cmd;  }  else{    var mv_cmd = mv + "somaticSV.vcf.gz " + inputs.output_basename + ".manta.somaticSV.vcf.gz" + mv + "somaticSV.vcf.gz.tbi " + inputs.output_basename + ".manta.somaticSV.vcf.gz.tbi" + mv + "candidateSmallIndels.vcf.gz " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz" + mv + "candidateSmallIndels.vcf.gz.tbi " + inputs.output_basename + ".manta.candidateSmallIndels.vcf.gz.tbi";    return "--tumorBam " + inputs.input_tumor_cram.path + " --normalBam " + inputs.input_normal_cram.path + std + mv_cmd;  }}</js> \
    """

}
