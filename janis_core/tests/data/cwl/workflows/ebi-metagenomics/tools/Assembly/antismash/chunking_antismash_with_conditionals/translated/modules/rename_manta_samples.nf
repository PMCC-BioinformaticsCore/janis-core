nextflow.enable.dsl=2

process RENAME_MANTA_SAMPLES {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest"
    publishDir "${params.outdir}/run_manta/rename_manta_samples"
    cpus "${params.run_manta.rename_manta_samples.cpus}"
    memory "${params.run_manta.rename_manta_samples.memory}"

    input:
    path input_vcf, stageAs: 'input_vcf'
    val input_normal_name
    val input_tumor_name

    output:
    tuple path("*.vcf.gz"), path("*.tbi"), emit: reheadered_vcf

    script:
    """
    echo \
    <js>inputs.input_normal_name) > sample_list.txt && echo $(inputs.input_tumor_name) >> sample_list.txt && bcftools reheader -s sample_list.txt $(inputs.input_vcf.path) > $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz")) && tabix $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz")</js> \
    """

}
