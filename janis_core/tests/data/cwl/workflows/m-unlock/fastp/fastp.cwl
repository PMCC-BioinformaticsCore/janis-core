#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

doc: |
  Modified from https://github.com/ambarishK/bio-cwl-tools/blob/release/fastp/fastp.cwl

requirements:
 - class: InlineJavascriptRequirement

hints:
  SoftwareRequirement:
    packages:
      fastp :
        version: ["0.23.2"]
        specs: ["https://anaconda.org/bioconda/fastp"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/fastp:0.23.2

inputs:
  identifier:
    type: string
    doc: Identifier for this dataset used in this workflow
    label: identifier used
  threads:
    type: int?
    default: 1
    inputBinding:
      prefix: --thread
  forward_reads:
    type: File
    inputBinding:
      prefix: --in1
  reverse_reads:
    type: File
    inputBinding:
      prefix: --in2
  merge_reads:
    type: boolean?
    default: false
    inputBinding:
      prefix: --merge
  qualified_phred_quality:
    type: int?
    default: 20
    inputBinding:
      prefix: --qualified_quality_phred
  unqualified_phred_quality:
    type: int?
    default: 20
    inputBinding:
      prefix: --unqualified_percent_limit
  min_length_required:
    type: int?
    default: 50
    inputBinding:
      prefix: --length_required
  force_polyg_tail_trimming:
    type: boolean?
    inputBinding:
      prefix: --trim_poly_g
  disable_trim_poly_g:
    type: boolean?
    default: true
    inputBinding:
      prefix: --disable_trim_poly_g
  base_correction:
    type: boolean?
    default: true
    inputBinding:
      prefix: --correction
  deduplicate:
    type: boolean?
    default: false
    inputBinding:
      prefix: --dedup

arguments:
  - prefix: --out1
    valueFrom: $(inputs.identifier)_fastp_1.fq.gz
  - |
    ${
      if (inputs.reverse_reads){
        return '--out2';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.reverse_reads){
        return inputs.identifier + "_fastp_2.fq.gz";
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.reverse_reads_path){
        return '--out2';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.reverse_reads_path){
        return inputs.identifier + "_fastp_2.fq.gz";
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.merge_reads){
        return '--merged_out';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.merge_reads){
        return inputs.identifier + "merged_fastp.fq.gz";
      } else {
        return '';
      }
    }
  
  - prefix: "-h"
    valueFrom: $(inputs.identifier)_fastp.html
  - prefix: "-j"
    valueFrom: $(inputs.identifier)_fastp.json


baseCommand: [fastp]

outputs:
  out_forward_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_fastp_1.fq.gz
  out_reverse_reads:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_fastp_2.fq.gz
  merged_reads:
    type: File?
    outputBinding:
      glob: $(inputs.identifier)_merged_fastp.fq.gz

  html_report:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_fastp.html
  json_report:
    type: File
    outputBinding:
      glob: $(inputs.identifier)_fastp.json

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
s:dateModified: "2022-02-22"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/