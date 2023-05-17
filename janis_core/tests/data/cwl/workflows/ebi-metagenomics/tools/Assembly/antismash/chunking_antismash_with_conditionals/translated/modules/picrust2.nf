nextflow.enable.dsl=2

process PICRUST2 {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/picrust2"

    input:
    path fasta, stageAs: 'fasta'
    path input_table, stageAs: 'input_table'

    output:
    path "{inputs.identifier}_PICRUSt2/COG_metagenome_out", emit: COG_metagenome_out
    path "{inputs.identifier}_PICRUSt2/EC_metagenome_out", emit: EC_metagenome_out
    path "{inputs.identifier}_PICRUSt2/EC_predicted.tsv.gz", emit: EC_predicted.tsv.gz
    path "{inputs.identifier}_PICRUSt2/KO_metagenome_out", emit: KO_metagenome_out
    path "{inputs.identifier}_PICRUSt2/KO_predicted.tsv.gz", emit: KO_predicted.tsv.gz
    path "{inputs.identifier}_PICRUSt2/PFAM_metagenome_out", emit: PFAM_metagenome_out
    path "{inputs.identifier}_PICRUSt2/PFAM_predicted.tsv.gz", emit: PFAM_predicted.tsv.gz
    path "{inputs.identifier}_PICRUSt2/TIGRFAM_metagenome_out", emit: TIGRFAM_metagenome_out
    path "{inputs.identifier}_PICRUSt2/TIGRFAM_predicted.tsv.gz", emit: TIGRFAM_predicted.tsv.gz
    path "{inputs.identifier}_PICRUSt2/intermediate", emit: intermediate
    path "{inputs.identifier}_PICRUSt2/marker_predicted_and_nsti.tsv.gz", emit: marker_predicted_and_nsti.tsv.gz
    path "{inputs.identifier}_PICRUSt2/out.tre", emit: out.tre
    path "{inputs.identifier}_PICRUSt2/pathways_out", emit: pathways_out
    path "{inputs.identifier}_picrust2.stderr.log", emit: stderr_out
    path "{inputs.identifier}_picrust2.stdout.log", emit: stdout_out

    script:
    def stratified = params.picrust2_stratified == false ? "" : "--stratified"
    def threads = params.threads ? params.threads : 2
    """
    picrust2_pipeline.py \
    -i ${input_table} \
    -s ${fasta} \
    -p ${threads} \
    ${stratified} \
    --in_traits COG,EC,KO,PFAM,TIGRFAM \
    --in_traits ${"COG,EC,KO,PFAM,TIGRFAM"} \
    -o "${params.sample}_PICRUSt2" \
    2> "${params.sample}_picrust2.stderr.log" \
    > "${params.sample}_picrust2.stdout.log" \
    """

}
