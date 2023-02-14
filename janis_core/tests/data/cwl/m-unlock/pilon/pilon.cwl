#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "Pilon"
doc: |
       "https://github.com/broadinstitute/pilon
        Pilon is a software tool which can be used to:
          Automatically improve draft assemblies
          Find variation among strains, including large event detection"

requirements:
 - class: InlineJavascriptRequirement 

hints:
  SoftwareRequirement:
    packages:
      pilon :
        version: ["1.24"]
        specs: ["https://anaconda.org/bioconda/pilon"]

  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/pilon:1.24

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used

  threads:
    type: int?
    doc: number of threads to use for computational processes
    label: number of threads
    inputBinding:
      prefix: --threads
    default: 2
  memory:
    type: int?
    doc: Memory usage in megabytes
    label: memory usage (MB)
    default: 8000

  assembly:
    type: File
    label: Assembly
    doc: Draft assembly fasta file
    inputBinding:
      prefix: --genome

  bam_file:
    type: File
    label: Bam file
    doc: Indexed sorted bam file with mapped reads to draft assembly
    inputBinding:
      prefix: --frags
    # secondaryFiles:
    #   - .bai
  fixlist:
    type: string?
    label: Fix List
    doc: | 
      A comma-separated list of categories of issues to try to fix:
        "snps": try to fix individual base errors;
        "indels": try to fix small indels;
        "gaps": try to fill gaps;
        "local": try to detect and fix local misassemblies;
        "all": all of the above (default);
        "bases": shorthand for "snps" and "indels" (for back compatibility);
    inputBinding:
      prefix: --fix

  vcf:
    label: vcf output
    type: boolean
    default: true
    inputBinding:
      prefix: --vcf

baseCommand: ["java"]

arguments:
  - "-jar"
  - "-Xmx$(inputs.memory)M"
  - "/venv/share/pilon-1.24-0/pilon.jar"
  - valueFrom: $(inputs.identifier)_pilon_polished
    prefix: "--output"


stdout: $(inputs.identifier)_pilon.log

outputs:
  pilon_polished_assembly:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_pilon_polished.fasta
  pilon_vcf:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_pilon_polished.vcf
  pilon_log:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_pilon.log

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
s:dateCreated: "2022-02-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
