#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: "spades genomic assembler"

doc: |
    Runs the spades assembler using a dataset file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: input_spades.json
        entry: |-
          [
            {
              orientation: "fr",
              type: "paired-end",
              right reads: $( inputs.forward_reads.map( function(x) {return  x.path} ) ),
              left reads: $( inputs.reverse_reads.map( function(x) {return  x.path} ) )
            }            
            ${
              var pacbio=""
                if (inputs.pacbio_reads != null) {
                 pacbio+=',{ type: "pacbio", single reads: ["' + inputs.pacbio_reads.map( function(x) {return  x.path} ).join('","') + '"] }' 
              }
              return pacbio;
            }
            ${
              var nanopore=""
                if (inputs.nanopore_reads != null) {
                 nanopore+=',{ type: "nanopore", single reads: ["' + inputs.nanopore_reads.map( function(x) {return  x.path} ).join('","') + '"] }'
                //  nanopore+=',{ type: "nanopore", single reads: ["' + inputs.nanopore_reads.join('","') + '"] }'
              }
              return nanopore;
            }
          ]

hints:
  SoftwareRequirement:
    packages:
      spades:
        version: ["3.15.5"]
        specs: ["https://anaconda.org/bioconda/spades"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/spades:3.15.5

# baseCommand: ["/unlock/infrastructure/binaries/SPAdes/SPAdes_v3.15.4/bin/spades.py", "--dataset", "input_spades.json"]
baseCommand: [spades.py, --dataset, input_spades.json]

arguments:
  - valueFrom: $(runtime.outdir)/output
    prefix: -o
  - valueFrom: $(inputs.memory / 1000)
    prefix: --memory

inputs:
  isolate:
    type: boolean?
    doc: this flag is highly recommended for high-coverage isolate and multi-cell data
    label: high-coverage mode
    inputBinding:
      prefix: --isolate
  metagenome:
    type: boolean?
    doc: this flag is required for metagenomic sample data
    label: metagenomics sample
    inputBinding:
      prefix: --meta
  biosyntheticSPAdes:
    type: boolean?
    doc: this flag is required for biosyntheticSPAdes mode
    label: biosynthetic spades mode
    inputBinding:
      prefix: --bio
  rna:
    type: boolean?
    doc: this flag is required for RNA-Seq data
    label: rnaseq data
    inputBinding:
      prefix: --rna
  plasmid:
    type: boolean?
    doc: runs plasmidSPAdes pipeline for plasmid detection
    label: plasmid spades run
    inputBinding:
      prefix: --plasmid
  only_assembler:
    type: boolean?
    doc: Runs only assembling (without read error correction)
    label: Only assembler
    inputBinding:
      prefix: --only-assembler

  IonTorrent:
    type: boolean?
    doc: this flag is required for IonTorrent data
    label: iontorrent data
    inputBinding:
      prefix: --iontorrent
  forward_reads:
    type: File[]
    doc: The file containing the forward reads
    label: Forward reads
  reverse_reads:
    type: File[]
    doc: The file containing the reverse reads
    label: Reverse reads
  pacbio_reads:
    type: File[]?
    doc: Fastq file with PacBio CLR reads
    label: PacBio CLR reads
  nanopore_reads:
    type: File[]?
    doc: Fastq file with Oxford NanoPore reads
    label: NanoPore reads
  threads:
    type: int
    doc: number of threads to use
    label: threads
    inputBinding:
      prefix: --threads
  memory:
    type: int?
    doc: Memory used in megabytes

outputs:
  # stdout: stdout
  # stderr: stderr
  # json:
  #   type: File
  #   outputBinding:
  #     glob: spades.json
  contigs:
    type: File
    outputBinding:
      glob: output/contigs.fasta

  scaffolds:
    type: File
    outputBinding:
      glob: output/scaffolds.fasta

  assembly_graph:
    type: File
    outputBinding:
      glob: output/assembly_graph.fastg

  contigs_assembly_paths:
    type: File
    outputBinding:
      glob: output/contigs.paths

  scaffolds_assembly_paths:
    type: File
    outputBinding:
      glob: output/scaffolds.paths

  contigs_before_rr:
    label: contigs before repeat resolution
    type: File
    outputBinding:
      glob: output/before_rr.fasta

  params:
    label: information about SPAdes parameters in this run
    type: File
    outputBinding:
      glob: output/params.txt

  log:
    label: SPAdes log
    type: File
    outputBinding:
      glob: output/spades.log

  internal_config:
    label: internal configuration file
    type: File
    outputBinding:
      glob: output/dataset.info

  internal_dataset:
    label: internal YAML data set file
    type: File
    outputBinding:
      glob: output/input_dataset.yaml

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
s:license: https://spdx.org/licenses/CC0-1.0.html 
s:copyrightHolder: "UNLOCK - Unlocking Microbial Potential"


$namespaces:
 s: http://schema.org/
