nextflow.enable.dsl=2

process EXTRACT_COORDS {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/other_ncrnas/extract_coords"
    cpus "${params.after_qc.other_ncrnas.extract_coords.cpus}"
    memory "${params.after_qc.other_ncrnas.extract_coords.memory}"

    input:
    path infernal_matches, stageAs: 'infernal_matches'
    val name

    output:
    path "*matched_seqs_with_coords*", emit: matched_seqs_with_coords

    script:
    def name = name ? name : 
    """
    awk_tool \
    -i ${infernal_matches} \
    -n ${name} \
    """

}
