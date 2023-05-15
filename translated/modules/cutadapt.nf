nextflow.enable.dsl=2

process CUTADAPT {
    
    container "quay.io/biocontainers/cutadapt:4.0--py37h8902056_0"
    publishDir "${params.outdir}/cutadapt"

    input:
    path unknown1
    path unknown2

    output:
    path "out1*", emit: out1
    path "out2*", emit: out2

    script:
    """
    cutadapt \
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
