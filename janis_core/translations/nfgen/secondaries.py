
from janis_core import ToolInput
from janis_core.types import (
    File,
    Array,
    DataType
)
from . import nfgen_utils


"""
def apply_secondary_file_format_to_filename()  EXISTS 
    - from janis_core.utils.secondary import apply_secondary_file_format_to_filename
    - applies janis ext conventions to create full filename (incl wildcards)
"""


def get_extensions(dtype: File) -> list[str]:
    """returns extension of each file for File types with secondaries"""
    primary_ext: str = ''
    secondary_exts: list[str] = []

    # primary extension
    if len(dtype.get_extensions()) > 0:
        primary_ext = dtype.get_extensions()[0]
    else:
        primary_ext = 'primary'
    
    # secondary extensions
    if dtype.secondary_files() is not None:
        secondary_exts = dtype.secondary_files()
    else:
        secondary_exts = []

    return sort_extensions(primary_ext, secondary_exts)

def get_names(dtype: File) -> list[str]:
    """returns name of each file for File types with secondaries"""
    exts = get_extensions(dtype)
    return remove_symbols(exts)
 
def sort_extensions(primary_ext: str, secondary_exts: list[str]) -> list[str]:
    out: list[str] = []
    out.append(primary_ext)
    secondary_exts = sorted(secondary_exts, key=lambda x: x.rsplit('.')[-1])
    out += secondary_exts
    return out

def remove_symbols(exts: list[str]) -> list[str]:
    return [x.rsplit('.')[-1] for x in exts]

def is_secondary_type(dtype: DataType) -> bool:
    basetype = nfgen_utils.get_base_type(dtype)
    if isinstance(basetype, File) and basetype.has_secondary_files():
        return True
    return False

def is_array_secondary_type(dtype: DataType) -> bool:
    if isinstance(dtype, Array) and is_secondary_type(dtype):
        return True
    return False




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