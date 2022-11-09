

"""
process GRIDSS {
  tag "gridss"
  container 'quay.io/biocontainers/gridss:2.9.3--0'

  input:
    tuple path(input), path(bai)
    tuple path(index_files), path(fai), path(reference_genome)
    path blacklist    

  output:
    path "*.gridss.raw.vcf.gz", emit: output

  script:
  def threadsArgument = params.threads_gridss ? "--threads $params.threads_gridss" : ""
  '''
  for i in $input
  do
    gridss --reference $reference_genome \\
            --output \$(echo \${i%%.*}).gridss.raw.vcf.gz \\
            --assembly "$params.assemblyFilename" \\
            $threadsArgument \\
            --jar "/usr/local/share/gridss-2.9.3-0/gridss.jar" \\
            --blacklist $blacklist \\
            \$i
  done
  '''
}

process SAMTOOLS_INDEX {
  tag "samtools index"
  container 'quay.io/biocontainers/samtools:1.5--2'

  input:
    path input

  output:
    tuple path(input), path("*.sorted*.bai")

  script:
  '''
  samtools index -b $input
  '''
}


process '6C_ASE_knownSNPs' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId"

  input:
    path genome from params.genome
    path index from genome_index_ch
    path dict from genome_dict_ch
    tuple val(sampleId), path(vcf), path(bam), path(bai) from grouped_vcf_bam_bai_ch
    
  output:
    path "ASE.tsv"

  script:
  '''
  echo "${bam.join('\n')}" > bam.list

  java -jar $GATK -R ${genome} \
                  -T ASEReadCounter \
                  -o ASE.tsv \
                  -I bam.list \
                  -sites ${vcf}
  '''
}



"""