cwlVersion: v1.2
class: CommandLineTool

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      bbmap:
        version: ["38.98"]
        specs: ["https://anaconda.org/bioconda/bbmap"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/bbmap:38.98
    
label: BBMap
doc: |
  Read filtering using BBMap against a (contamination) reference genome

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  forward_reads:
    type: File
    inputBinding:
      position: 1
      prefix: 'in='
      separate: false
  reverse_reads:
    type: File?
    inputBinding:
      position: 2
      prefix: 'in2='
      separate: false
  reference:
    type: File
    inputBinding:
      position: 3
      prefix: 'ref='
      separate: false
  output_mapped:
    type: boolean?
    default: false

  memory:
    type: int?
    doc: maximum memory usage in megabytes
    label: memory usage (mb)
    default: 8000
  threads:
    type: int?
    doc: number of threads to use for computational processes
    label: number of threads
    inputBinding:
      prefix: 'threads='
      separate: false
    default: 2

stderr: "$(inputs.identifier)_BBMap_log.txt"

outputs:
  log:
    label: BBMap log output
    type: File
    outputBinding:
      glob: "$(inputs.identifier)_BBMap_log.txt"
  stats:
    label: Mapping statistics
    type: File
    outputBinding:
      glob: "$(inputs.identifier)_BBMap_stats.txt"
  covstats:
    label: Coverage per contig
    type: File
    outputBinding:
      glob: "$(inputs.identifier)_BBMap_covstats.txt"
  out_forward_reads:
    type: File
    outputBinding:
      glob: "$(inputs.identifier)_filtered_1.fq.gz"
  out_reverse_reads:
    type: File
    outputBinding:
      glob: "$(inputs.identifier)_filtered_2.fq.gz"


baseCommand: [bbmap.sh]

arguments:
  - "-Xmx$(inputs.memory)M"
  - "printunmappedcount"
  - "overwrite=true"
  - "bloom=t"
  - "statsfile=$(inputs.identifier)_BBMap_stats.txt"
  - "covstats=$(inputs.identifier)_BBMap_covstats.txt"
  - |
    ${
      if (inputs.output_mapped){
        return 'outm1='+inputs.identifier+'_filtered_1.fq.gz \
                outm2='+inputs.identifier+'_filtered_2.fq.gz';
      } else {
        return 'outu1='+inputs.identifier+'_filtered_1.fq.gz \
                outu2='+inputs.identifier+'_filtered_2.fq.gz';
      }
    }
  # - "fast"
  # - "minratio=0.9"
  # - "maxindel=3"
  # - "bwr=0.16"
  # - "bw=12"
  # - "minhits=2"
  # - "qtrim=r"
  # - "trimq=10"
  # - "untrim"
  # - "idtag"
  # - "kfilter=25"
  # - "maxsites=1"
  # - "k=14"
  # - "nodisk=t"
  # - "out=$(inputs.identifier)_BBMap.sam"
  # - "rpkm=$(inputs.identifier).rpkm"

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-8172-8981
    s:email: mailto:jasper.koehorst@wur.nl
    s:name: Jasper Koehorst
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-9524-5964
    s:email: mailto:bart.nijsse@wur.nl
    s:name: Bart Nijsse

s:citation: https://m-unlock.nl
s:codeRepository: https://gitlab.com/m-unlock/cwl
s:dateCreated: "2020-00-00"
s:dateModified: "2022-04-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/