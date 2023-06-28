nextflow.enable.dsl=2

process CMSEARCH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/cmsearch"
    cpus "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch.cpus}"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch.memory}"

    input:
    path query_sequences, stageAs: 'query_sequences'
    val covariance_model_database

    output:
    path "<js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch_matches.tbl";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch_matches.tbl";  }  return name;}</js>", emit: matches
    path "<js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch.out";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch.out";  }  return name;}</js>", emit: programOutput

    script:
    def cpu = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cpu ? "--cpu ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cpu}" : ""
    def cut_ga = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cut_ga ? "--cut_ga" : ""
    def omit_alignment_section = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_omit_alignment_section == false ? "" : "--noali"
    def only_hmm = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_only_hmm ? "--hmmonly" : ""
    """
    cmsearch \
    ${cpu} \
    -Z ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_search_space_size} \
    ${cut_ga} \
    ${only_hmm} \
    ${omit_alignment_section} \
    --tblout <js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch_matches.tbl";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch_matches.tbl";  }  return name;}</js> \
    -o <js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch.out";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch.out";  }  return name;}</js> \
    ${covariance_model_database} \
    ${query_sequences} \
    > /dev/null \
    2> /dev/null \
    """

}


process CMSEARCH {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/after_qc/rna_prediction/find_ribosomal_ncrnas/cmsearch_raw_data/cmsearch"
    cpus "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch.cpus}"
    memory "${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch.memory}"

    input:
    path query_sequences, stageAs: 'query_sequences'
    val covariance_model_database

    output:
    path "<js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch_matches.tbl";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch_matches.tbl";  }  return name;}</js>", emit: matches
    path "<js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch.out";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch.out";  }  return name;}</js>", emit: programOutput

    script:
    def cpu = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cpu ? "--cpu ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cpu}" : ""
    def cut_ga = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_cut_ga ? "--cut_ga" : ""
    def omit_alignment_section = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_omit_alignment_section == false ? "" : "--noali"
    def only_hmm = params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_only_hmm ? "--hmmonly" : ""
    """
    cmsearch \
    ${cpu} \
    -Z ${params.after_qc.rna_prediction.find_ribosomal_nc_rnas.cmsearch_raw_data.cmsearch_search_space_size} \
    ${cut_ga} \
    ${only_hmm} \
    ${omit_alignment_section} \
    --tblout <js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch_matches.tbl";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch_matches.tbl";  }  return name;}</js> \
    -o <js>${  var name = "";  if (typeof inputs.covariance_model_database == "string") {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.split("/").slice(-1)[0] +      ".cmsearch.out";  } else {    name =      inputs.query_sequences.basename +      "." +      inputs.covariance_model_database.nameroot +      ".cmsearch.out";  }  return name;}</js> \
    ${covariance_model_database} \
    ${query_sequences} \
    > /dev/null \
    2> /dev/null \
    """

}
