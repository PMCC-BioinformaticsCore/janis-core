nextflow.enable.dsl=2

include { RETURN_ANTISMASH_IN_FOLDER } from '../modules/return_antismash_in_folder'
include { TOUCH_NO_ANTISMASH_FLAG } from '../modules/touch_no_antismash_flag'

workflow NO_ANTISMASH_SUBWF {

    take:
    ch_filtered_fasta
    ch_final_folder_name

    main:
    RETURN_ANTISMASH_IN_FOLDER(
        TOUCH_NO_ANTISMASH_FLAG.out.created_file.toList(),
        ch_final_folder_name
    )

    TOUCH_NO_ANTISMASH_FLAG()

    emit:
    antismash_result_folder = RETURN_ANTISMASH_IN_FOLDER.out.out

}
