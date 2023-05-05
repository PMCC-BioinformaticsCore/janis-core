nextflow.enable.dsl=2

process PROKKA {
    debug true
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
    --cpus 8 \
    --evalue 1e-06 \
    --gcode 11 \
    --genus "Escherichia" \
    --gffver 3 \
    --increment 10 \
    --kingdom "Bacteria" \
    --locustag "PROKKA" \
    --mincontig 200 \
    --outdir "outdir" \
    --prefix "prokka" \
    --species "Coli" \
    --strain "C-1" \
    ${input_file} \
    """

}
