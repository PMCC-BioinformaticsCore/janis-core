nextflow.enable.dsl=2

process GATK_SELECTVARIANTS_MANTA {
    debug true
    container "pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0"
    publishDir "${params.outdir}/run_manta/gatk_selectvariants_manta"
    cpus "${params.run_manta.gatk_selectvariants_manta.cpus}"
    memory "${params.run_manta.gatk_selectvariants_manta.memory}"

    input:
    tuple path(primary), path(tbi)
    val output_basename
    val mode

    output:
    tuple path("*.PASS.vcf.gz"), path("*.tbi"), emit: pass_vcf

    script:
    """
    /bin/bash -c \
    set -eo pipefail
    ${
      var run_mode = inputs.mode;
      if (run_mode == 'grep' || run_mode == 'gatk'){
        var in_vcf = inputs.input_vcf.path;
        var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
        var cmd = '/gatk SelectVariants --java-options "-Xmx8000m" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
        if (run_mode == 'grep'){
          cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
        }
        return cmd;
      }
      else{
        throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
      }
    } \
    """

}
