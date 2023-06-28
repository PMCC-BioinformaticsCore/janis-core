#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'BCFTools: Annotate'
doc: "------------------------------------\n\nAdd or remove annotations."

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biocontainers/bcftools:v1.5_cv2

inputs:
- id: vcf
  label: vcf
  type: File
  inputBinding:
    position: 10
- id: outputFilename
  label: outputFilename
  doc: '[-o] see Common Options'
  type:
  - string
  - 'null'
  default: generated.vcf
  inputBinding:
    prefix: --output
- id: annotations
  label: annotations
  doc: |-
    [-a] Bgzip-compressed and tabix-indexed file with annotations. The file can be VCF, BED, or a tab-delimited file with mandatory columns CHROM, POS (or, alternatively, FROM and TO), optional columns REF and ALT, and arbitrary number of annotation columns. BED files are expected to have the ".bed" or ".bed.gz" suffix (case-insensitive), otherwise a tab-delimited file is assumed. Note that in case of tab-delimited file, the coordinates POS, FROM and TO are one-based and inclusive. When REF and ALT are present, only matching VCF records will be annotated. When multiple ALT alleles are present in the annotation file (given as comma-separated list of alleles), at least one must match one of the alleles in the corresponding VCF record. Similarly, at least one alternate allele from a multi-allelic VCF record must be present in the annotation file. Missing values can be added by providing "." in place of actual value. Note that flag types, such as "INFO/FLAG", can be annotated by including a field with the value "1" to set the flag, "0" to remove it, or "." to keep existing flags. See also -c, --columns and -h, --header-lines.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --annotations
- id: collapse
  label: collapse
  doc: |-
    (snps|indels|both|all|some|none) Controls how to match records from the annotation file to the target VCF. Effective only when -a is a VCF or BCF. See Common Options for more.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --collapse
- id: columns
  label: columns
  doc: |-
    [-c] Comma-separated list of columns or tags to carry over from the annotation file (see also -a, --annotations). If the annotation file is not a VCF/BCF, list describes the columns of the annotation file and must include CHROM, POS (or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which should be ignored can be indicated by "-". If the annotation file is a VCF/BCF, only the edited columns/tags must be present and their order does not matter. The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags can be written both as "INFO/TAG" or simply "TAG", and FORMAT tags can be written as "FORMAT/TAG" or "FMT/TAG". The imported VCF annotations can be renamed as "DST_TAG:=SRC_TAG" or "FMT/DST_TAG:=FMT/SRC_TAG". To carry over all INFO annotations, use "INFO". To add all INFO annotations except "TAG", use "^INFO/TAG". By default, existing values are replaced. To add annotations without overwriting existing values (that is, to add missing tags or add values to existing tags with missing values), use "+TAG" instead of "TAG". To append to existing values (rather than replacing or leaving untouched), use "=TAG" (instead of "TAG" or "+TAG"). To replace only existing values without modifying missing annotations, use "-TAG". If the annotation file is not a VCF/BCF, all new annotations must be defined via -h, --header-lines.
  type:
  - type: array
    items: string
  - 'null'
  inputBinding:
    prefix: --columns
- id: exclude
  label: exclude
  doc: |-
    [-e] exclude sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --exclude
- id: headerLines
  label: headerLines
  doc: |-
    [-h] Lines to append to the VCF header, see also -c, --columns and -a, --annotations.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --header-lines
- id: setId
  label: setId
  doc: |-
    [-I] assign ID on the fly. The format is the same as in the query command (see below). By default all existing IDs are replaced. If the format string is preceded by "+", only missing IDs will be set. For example, one can use # bcftools annotate --set-id +' % CHROM\_ % POS\_ % REF\_ % FIRST_ALT' file.vcf
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --set-id
- id: include
  label: include
  doc: |-
    [-i] include only sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --include
- id: keepSites
  label: keepSites
  doc: keep sites wich do not pass -i and -e expressions instead of discarding them(
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --keep-sites
- id: markSites
  label: markSites
  doc: |-
    [-m] (+|-)annotate sites which are present ("+") or absent ("-") in the -a file with a new INFO/TAG flag
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --mark-sites
- id: outputType
  label: outputType
  doc: '[-O] (b|u|z|v) see Common Options'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --output-type
- id: regions
  label: regions
  doc: ([-r] chr|chr:pos|chr:from-to|chr:from-[,…]) see Common Options
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --regions
- id: regionsFile
  label: regionsFile
  doc: '[-R] see Common Options'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --regions-file
- id: renameChrs
  label: renameChrs
  doc: |-
    rename chromosomes according to the map in file, with "old_name new_name\n" pairs separated by whitespaces, each on a separate line.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --rename-chrs
- id: samples
  label: samples
  doc: '[-s] subset of samples to annotate, see also Common Options'
  type:
  - type: array
    items: File
  - 'null'
  inputBinding:
    prefix: --samples
- id: samplesFile
  label: samplesFile
  doc: |-
    [-S] subset of samples to annotate. If the samples are named differently in the target VCF and the -a, --annotations VCF, the name mapping can be given as "src_name dst_name\n", separated by whitespaces, each pair on a separate line.
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --samples-file
- id: threads
  label: threads
  doc: see Common Options
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --threads
- id: remove
  label: remove
  doc: |-
    [-x] List of annotations to remove. Use "FILTER" to remove all filters or "FILTER/SomeFilter" to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and "FORMAT" to remove all FORMAT tags except GT. To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be abbreviated to "INF" and "FORMAT" to "FMT".
  type:
  - type: array
    items: string
  - 'null'
  inputBinding:
    prefix: --remove

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- bcftools
- annotate
arguments: []
id: bcftoolsAnnotate
