nextflow.enable.dsl=2

include { GATK_SELECTVARIANTS_MANTA } from '../modules/gatk_selectvariants_manta'
include { MANTA } from '../modules/manta'
include { RENAME_MANTA_SAMPLES } from '../modules/rename_manta_samples'

workflow RUN_MANTA {

    take:
    ch_hg38_strelka_bed
    ch_indexed_reference_fasta
    ch_input_normal_aligned
    ch_input_normal_name
    ch_input_tumor_aligned
    ch_input_tumor_name
    ch_output_basename
    ch_reference_dict
    ch_vep_cache
    ch_manta_cores
    ch_manta_memory
    ch_select_vars_mode

    main:
    GATK_SELECTVARIANTS_MANTA(
        RENAME_MANTA_SAMPLES.out.reheadered_vcf,
        ch_output_basename,
        ch_select_vars_mode
    )

    MANTA(
        ch_hg38_strelka_bed,
        ch_input_normal_aligned,
        ch_input_tumor_aligned,
        ch_indexed_reference_fasta,
        ch_output_basename,
        ch_manta_cores,
        ch_manta_memory
    )

    RENAME_MANTA_SAMPLES(
        MANTA.out.output_sv.map{ tuple -> tuple[0] },
        ch_input_normal_name,
        ch_input_tumor_name
    )

    emit:
    manta_pass_vcf = GATK_SELECTVARIANTS_MANTA.out.pass_vcf
    manta_prepass_vcf = RENAME_MANTA_SAMPLES.out.reheadered_vcf
    manta_small_indels = MANTA.out.small_indels

}
