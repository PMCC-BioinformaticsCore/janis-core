nextflow.enable.dsl=2

process PILON {
    debug true
    container "ubuntu:latest"
    publishDir "${params.outdir}/pilon"

    input:
    path assembly, stageAs: 'assembly'
    path bam_file, stageAs: 'bam_file'

    output:
    path "{inputs.identifier}_pilon.log", emit: pilon_log
    path "{inputs.identifier}_pilon_polished.fasta", emit: pilon_polished_assembly
    path "{inputs.identifier}_pilon_polished.vcf", emit: pilon_vcf

    script:
    def fixlist = params.fixlist ? "--fix ${params.fixlist}" : ""
    def threads = params.threads ? params.threads : 2
    """
    java \
    --frags ${bam_file} \
    --genome ${assembly} \
    ${fixlist} \
    --threads ${threads} \
    --vcf \
    --output "${params.identifier}_pilon_polished" \
    -jar \
    "-Xmx${params.memory}M" \
    /venv/share/pilon-1.24-0/pilon.jar \
    > "${params.identifier}_pilon.log" \
    """

}
