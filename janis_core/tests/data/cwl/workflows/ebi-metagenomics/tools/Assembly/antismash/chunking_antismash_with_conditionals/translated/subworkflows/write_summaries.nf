nextflow.enable.dsl=2

include { FUNCTIONAL_STATS } from '../modules/functional_stats'
include { WRITE_SUMMARIES } from '../modules/write_summaries'

workflow WRITE_SUMMARIES {

    take:
    ch_cds
    ch_hmmscan_annotation
    ch_interproscan_annotation
    ch_ko_file
    ch_pfam_annotation
    ch_rna

    main:
    FUNCTIONAL_STATS(
        ch_cds,
        ch_rna,
        ch_hmmscan_annotation,
        ch_interproscan_annotation,
        ch_pfam_annotation,
        ch_ko_file
    )

    WRITE_SUMMARIES(
        FUNCTIONAL_STATS.out.ips_yaml,
        FUNCTIONAL_STATS.out.ko_yaml,
        FUNCTIONAL_STATS.out.pfam_yaml,
        FUNCTIONAL_STATS.out.antismash_yaml,
        ch_cds,
        ch_cds,
        ch_cds,
        ch_cds
    )

    emit:
    stats = FUNCTIONAL_STATS.out.stats
    summary_antismash = WRITE_SUMMARIES.out.summary_antismash
    summary_ips = WRITE_SUMMARIES.out.summary_ips
    summary_ko = WRITE_SUMMARIES.out.summary_ko
    summary_pfam = WRITE_SUMMARIES.out.summary_pfam

}
