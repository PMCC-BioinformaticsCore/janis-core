nextflow.enable.dsl=2

process PROKKA {
    
    container "quay.io/biocontainers/prokka:1.14.5--pl526_1"
    publishDir "${params.outdir}/prokka"

    input:
    path input_file

    output:
    path "outdir/prokka.err", emit: outErr
    path "outdir/prokka.faa", emit: outFaa
    path "outdir/prokka.ffn", emit: outFfn
    path "outdir/prokka.fna", emit: outFna
    path "outdir/prokka.fsa", emit: outFsa
    path "outdir/prokka.gbk", emit: outGbk
    path "outdir/prokka.gff", emit: outGff
    path "outdir/prokka.log", emit: outLog
    path "outdir/prokka.sqn", emit: outSqn
    path "outdir/prokka.tbl", emit: outTbl
    path "outdir/prokka.txt", emit: outTxt

    script:
    """
    prokka \
    ${input_file} \
    """

}
