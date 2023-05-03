nextflow.enable.dsl=2

process CUTADAPT {
    debug true
    container "quay.io/biocontainers/cutadapt:4.0--py37h8902056_0"
    publishDir "${params.outdir}/cutadapt"

    input:
    path output_file
    path unknown1
    path unknown2

    output:
    path "out1*", emit: out11
    path "out1*", emit: out12
    path "out2*", emit: out21
    path "out2*", emit: out22
    path "out1.fq*", emit: outForward
    path None, emit: outInfoFile
    path "stats.json", emit: outJsonStats
    path "report.txt", emit: outReport
    path "rest_output*", emit: outRestOutput1
    path "rest_output*", emit: outRestOutput2
    path "out2.fq*", emit: outReverse
    path "split/.+\.fastq.*", emit: outSplitOutput
    path "too_long_output*", emit: outTooLongOutput1
    path "too_long_output*", emit: outTooLongOutput2
    path "too_long_paired_output*", emit: outTooLongPairedOutput1
    path "too_long_paired_output*", emit: outTooLongPairedOutput2
    path "too_short_output*", emit: outTooShortOutput1
    path "too_short_output*", emit: outTooShortOutput2
    path "too_short_paired_output*", emit: outTooShortPairedOutput1
    path "too_short_paired_output*", emit: outTooShortPairedOutput2
    path "untrimmed_output*", emit: outUntrimmedOutput1
    path "untrimmed_output*", emit: outUntrimmedOutput2
    path "untrimmed_paired_output*", emit: outUntrimmedPairedOutput1
    path "untrimmed_paired_output*", emit: outUntrimmedPairedOutput2
    path "wild_output*", emit: outWildOutput1
    path "wild_output*", emit: outWildOutput2

    script:
    def info_file = null
    """
    cutadapt \
    --output=${output_file} \
    --action="trim" \
    --error-rate 0.1 \
    --json "stats.json" \
    --length 0 \
    --nextseq-trim 0 \
    --overlap 3 \
    --pair-filter "any" \
    --quality-cutoff "0" \
    --times 1 \
    --zero-cap "" \
    -U 0 \
    -j=4 \
    -u 0 \
    """

}
