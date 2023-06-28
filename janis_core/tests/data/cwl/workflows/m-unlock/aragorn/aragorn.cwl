#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: "tRNA prediction using Aragorn"

requirements:
  ResourceRequirement:
    ramMin: 5000
    coresMin: 4

hints:
  SoftwareRequirement:
    packages:
      aragorn:
        version: ["1.2.41"]
        specs: ["https://anaconda.org/bioconda/aragorn"]
  DockerRequirement:
    dockerPull: docker-registry.wur.nl/m-unlock/docker/aragorn:1.2.41

# baseCommand: [ /unlock/infrastructure/binaries/aragorn/aragorn_v1.2.38/aragorn, -m, -t, ....]
baseCommand: [aragorn, -m, -t]

# String command = aragorn + " -w -fasta " + fastaFile + " -o " + fastaFile + "_aragorn.txt";

arguments:
  - valueFrom: "-w"
  - valueFrom: $(inputs.input_fasta.basename).aragorn.fasta
    prefix: "-o"
  - valueFrom: $(inputs.input_fasta.basename).prodigal.ffn
    prefix: "-d"
  - valueFrom: $(inputs.input_fasta.basename).prodigal.faa
    prefix: "-a"

inputs:
  input_fasta:
    type: File
    inputBinding:
      separate: true
      prefix: "-i"
  meta:
    type: boolean?
    doc: Input is a meta-genome
    inputBinding:
      prefix: -p meta
  single:
    type: boolean?
    doc: Input is an isolate genome
    inputBinding:
      prefix: -p single

stdout: stdout.txt
stderr: stderr.txt

outputs:
  stdout: stdout
  stderr: stderr

  predicted_proteins_out:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.basename).prodigal
  predicted_proteins_ffn:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.basename).prodigal.ffn
  predicted_proteins_faa:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.basename).prodigal.faa

's:author': 'Ekaterina Sakharova'
's:copyrightHolder': EMBL - European Bioinformatics Institute
's:license': "https://www.apache.org/licenses/LICENSE-2.0"

$namespaces:
 s: http://schema.org/
