#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
requirements:
   - class: StepInputExpressionRequirement
   - class: InlineJavascriptRequirement
   - class: MultipleInputFeatureRequirement

label: Genome conversion and annotation
doc: Workflow for genome annotation from EMBL format

inputs:
  embl:
    type: File
    doc: Genome sequence in EMBL format
    label: EMBL input file
  identifier:
    type: string
    doc: Identifier of the sample being converted
    label: Sample name
  codon:
    type: int
    doc: Codon table used for gene prediction
    label: Codon table
  threads:
    type: int?
    doc: number of threads to use for computational processes
    label: number of threads
    default: 2
  interpro:
    type: Directory
    doc: Path to the interproscan application directory 
    label: InterProScan path
  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.

steps:
  conversion:
    run: ../sapp/conversion.cwl
    in:
      identifier: identifier
      embl: embl
      codon: codon
    out: [output]
############################
  # infernal:
  #   run: ../sapp/infernal.cwl
  #   in:
  #     identifier: identifier
  #     input: conversion/output
  #     threads: threads
  #   out: [output]
############################
  prodigal:
    run: ../sapp/prodigal.cwl
    in:
      identifier: identifier
      input: conversion/output
    out: [output]
############################
  kofamscan:
    run: ../sapp/kofamscan.cwl
    in:
      identifier: identifier
      input: prodigal/output
      threads: threads
    out: [output]
############################
  interproscan:
    run: ../sapp/interproscan.cwl
    in:
      identifier: identifier
      input: kofamscan/output
      cpu: threads
      interpro: interpro
    out: [output]
############################
  gzip:
    run: ../bash/compress.cwl
    in:
      infile: interproscan/output
    out: [output]

outputs:
  output:
    type: File
    outputSource: gzip/output
  
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
s:dateCreated: "2021-00-00"
s:dateModified: "2022-05-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"

$namespaces:
  s: https://schema.org/
