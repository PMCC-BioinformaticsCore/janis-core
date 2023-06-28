nextflow.enable.dsl=2

process MEMOTE_REPORT_SNAPSHOT {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/memote_report_snapshot"

    input:
    path gem, stageAs: 'gem'

    output:
    path "{inputs.GEM.basename}_MEMOTE.html", optional: true, emit: report_html
    path "{inputs.GEM.basename}_MEMOTE.json.gz", optional: true, emit: run_json

    script:
    def report_snapshot = params.memote_report_snapshot_report_snapshot == false ? "" : "None"
    def skip_test_find_metabolites_consumed_with_closed_bounds = params.memote_report_snapshot_skip_test_find_metabolites_consumed_with_closed_bounds == false ? "" : "--skip test_find_metabolites_consumed_with_closed_bounds"
    def skip_test_find_metabolites_not_consumed_with_open_bounds = params.memote_report_snapshot_skip_test_find_metabolites_not_consumed_with_open_bounds == false ? "" : "--skip test_find_metabolites_not_consumed_with_open_bounds"
    def skip_test_find_metabolites_not_produced_with_open_bounds = params.memote_report_snapshot_skip_test_find_metabolites_not_produced_with_open_bounds == false ? "" : "--skip test_find_metabolites_not_produced_with_open_bounds"
    def skip_test_find_metabolites_produced_with_closed_bounds = params.memote_report_snapshot_skip_test_find_metabolites_produced_with_closed_bounds == false ? "" : "--skip test_find_metabolites_produced_with_closed_bounds"
    """
    bash script.sh \
    ${report_snapshot} \
    <js>${  if (inputs.run){    return "run --filename " + inputs.GEM.basename + "_MEMOTE.json.gz";  } else {    return '';  }}</js> \
    <js>${  if (inputs.report_snapshot){    return "report snapshot --filename " + inputs.GEM.basename + "_MEMOTE.html";  } else {    return '';  }}</js> \
    ${skip_test_find_metabolites_produced_with_closed_bounds} \
    ${skip_test_find_metabolites_consumed_with_closed_bounds} \
    ${skip_test_find_metabolites_not_produced_with_open_bounds} \
    ${skip_test_find_metabolites_not_consumed_with_open_bounds} \
    ${gem} \
    """

}
