#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
requirements:
   - class: StepInputExpressionRequirement
   - class: InlineJavascriptRequirement
   - class: MultipleInputFeatureRequirement

inputs:
  forward_primer: 
    type: string
    doc: Forward primer used
    label: Forward primer
  reverse_primer:
    type: string?
    doc: Reverse primer used
    label: Reverse primer
    # default: "CGAC[AG][AG]CCATGCA[ACGT]CACCT"
  reference_db:
    type: string?
    doc: Reference database used in FASTA format
    label: Reference database
    # default: "/tempZone/References/Databases/Silva/SILVA_132_SSURef_tax_silva.fasta.gz"
  rev_read_len: 
    type: int?
    doc: Read length of the reverse read
    label: Reverse read length
  for_read_len: 
    type: int
    doc: Read length of the reverse read
    label: Reverse read length
  sample:
    type: string
    doc: Name of the sample being analysed
    label: Sample name
  mock3:
    type: string?
    doc: Mock3 reference selection
    label: Mock3 reference
  mock4:
    type: string?
    doc: Mock4 reference selection
    label: Mock4 reference
  fragment:
    type: string
    doc: Subfragment that is being analysed (e.g. V1-V3 or V5-region)
    label: Subfragment name
  metadata:
    type: File?
    doc: UNLOCK assay metadata file
    label: Metadata file

  destination:
    type: string?
    label: Output Destination
    doc: Optional Output destination used for cwl-prov reporting.
    
steps:
############################
  ngtax:
    run: ../ngtax/ngtax.cwl
    in:
      forward_primer: forward_primer
      reverse_primer: reverse_primer
      reference_db: reference_db
      mock3: mock3
      mock4: mock4
      rev_read_len: rev_read_len
      for_read_len: for_read_len
      sample: sample
      fragment: fragment
    out: [biom, turtle]
#############################
  ngtax_to_tsv-fasta:
    run: ../ngtax/ngtax_to_tsv-fasta.cwl
    in:
        input: ngtax/turtle
        identifier: sample
        fragment: fragment
        metadata: metadata
    out:
      [picrust_fasta, picrust_tsv, physeq_asv, physeq_seq, physeq_tax, physeq_met]
############################
  ngtax_files_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files:
        source: [ngtax/biom, ngtax/turtle]
      destination: 
        valueFrom: $("1_Classification")
    out:
      [results]
############################
  phyloseq_files_to_folder:
    run: ../expressions/files_to_folder.cwl
    in:
      files: 
        source: [ngtax_to_tsv-fasta/physeq_asv, ngtax_to_tsv-fasta/physeq_seq, ngtax_to_tsv-fasta/physeq_tax, ngtax_to_tsv-fasta/physeq_met]
        linkMerge: merge_flattened
      destination:
        valueFrom: $("2_PHYLOSEQ")
    out:
      [results]

outputs:
  files_to_folder_ngtax:
    type: Directory
    outputSource: ngtax_files_to_folder/results
  files_to_folder_phyloseq:
    type: Directory
    outputSource: phyloseq_files_to_folder/results


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
s:dateModified: "2022-05-00"
s:license: https://spdx.org/licenses/Apache-2.0 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
  s: https://schema.org/
