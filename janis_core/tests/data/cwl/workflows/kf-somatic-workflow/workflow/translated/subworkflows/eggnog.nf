nextflow.enable.dsl=2

include { EGGNOG_ANNOTATION } from '../modules/eggnog_annotation'
include { EGGNOG_HOMOLOGY_SEARCHES } from '../modules/eggnog_homology_searches'
include { UNITE_SEED_ORTHOLOGS } from '../modules/unite_seed_orthologs'

workflow EGGNOG {

    take:
    ch_cpu
    ch_data_dir
    ch_db
    ch_db_diamond
    ch_fasta_file
    ch_file_acc

    main:
    EGGNOG_ANNOTATION(
        UNITE_SEED_ORTHOLOGS.out.result,
        ch_cpu,
        ch_data_dir,
        ch_file_acc
    )

    EGGNOG_HOMOLOGY_SEARCHES(
        ch_fasta_file.flatten().first(),
        ch_cpu,
        ch_data_dir,
        ch_db,
        ch_db_diamond,
        ch_file_acc
    )

    UNITE_SEED_ORTHOLOGS(
        EGGNOG_HOMOLOGY_SEARCHES.out.output_orthologs.toList(),
        ch_file_acc
    )

    emit:
    annotations = EGGNOG_ANNOTATION.out.output_annotations
    orthologs = UNITE_SEED_ORTHOLOGS.out.result

}
