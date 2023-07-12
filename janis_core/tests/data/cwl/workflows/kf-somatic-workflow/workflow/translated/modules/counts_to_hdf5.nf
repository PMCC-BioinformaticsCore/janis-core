nextflow.enable.dsl=2

process COUNTS_TO_HDF5 {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/counts_to_hdf5"

    input:
    path biom, stageAs: 'biom'

    output:
    path "<js>${ var ext = "";if (inputs.json) { ext = "_json.biom"; }if (inputs.hdf5) { ext = "_hdf5.biom"; }if (inputs.tsv) { ext = "_tsv.biom"; }var pre = inputs.biom.nameroot.split('.');pre.pop()return pre.join('.') + ext; }</js>", emit: result

    script:
    def biom = biom ? "--input-fp ${biom}" : ""
    def hdf5 = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_hdf5 == false ? "" : "--to-hdf5"
    def table_type = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type ? "--table-type ${params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type}" : ""
    """
    biom-convert.sh \
    ${biom} \
    ${table_type} \
    ${hdf5} \
    --output-fp <js>${ var ext = "";   if (inputs.json) { ext = "_json.biom"; }   if (inputs.hdf5) { ext = "_hdf5.biom"; }   if (inputs.tsv) { ext = "_tsv.biom"; }   var pre = inputs.biom.nameroot.split('.');   pre.pop()   return pre.join('.') + ext; }</js> \
    """

}


process COUNTS_TO_HDF5 {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/counts_to_hdf5"

    input:
    path biom, stageAs: 'biom'

    output:
    path "<js>${ var ext = "";if (inputs.json) { ext = "_json.biom"; }if (inputs.hdf5) { ext = "_hdf5.biom"; }if (inputs.tsv) { ext = "_tsv.biom"; }var pre = inputs.biom.nameroot.split('.');pre.pop()return pre.join('.') + ext; }</js>", emit: result

    script:
    def biom = biom ? "--input-fp ${biom}" : ""
    def hdf5 = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_hdf5 == false ? "" : "--to-hdf5"
    def table_type = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type ? "--table-type ${params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type}" : ""
    """
    biom-convert.sh \
    ${biom} \
    ${table_type} \
    ${hdf5} \
    --output-fp <js>${ var ext = "";   if (inputs.json) { ext = "_json.biom"; }   if (inputs.hdf5) { ext = "_hdf5.biom"; }   if (inputs.tsv) { ext = "_tsv.biom"; }   var pre = inputs.biom.nameroot.split('.');   pre.pop()   return pre.join('.') + ext; }</js> \
    """

}


process COUNTS_TO_HDF5 {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/classify_ssus/counts_to_hdf5"

    input:
    path biom, stageAs: 'biom'

    output:
    path "<js>${ var ext = "";if (inputs.json) { ext = "_json.biom"; }if (inputs.hdf5) { ext = "_hdf5.biom"; }if (inputs.tsv) { ext = "_tsv.biom"; }var pre = inputs.biom.nameroot.split('.');pre.pop()return pre.join('.') + ext; }</js>", emit: result

    script:
    def biom = biom ? "--input-fp ${biom}" : ""
    def hdf5 = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_hdf5 == false ? "" : "--to-hdf5"
    def table_type = params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type ? "--table-type ${params.after_qc.rna_prediction.classify_ssus.counts_to_hdf5_table_type}" : ""
    """
    biom-convert.sh \
    ${biom} \
    ${table_type} \
    ${hdf5} \
    --output-fp <js>${ var ext = "";   if (inputs.json) { ext = "_json.biom"; }   if (inputs.hdf5) { ext = "_hdf5.biom"; }   if (inputs.tsv) { ext = "_tsv.biom"; }   var pre = inputs.biom.nameroot.split('.');   pre.pop()   return pre.join('.') + ext; }</js> \
    """

}
