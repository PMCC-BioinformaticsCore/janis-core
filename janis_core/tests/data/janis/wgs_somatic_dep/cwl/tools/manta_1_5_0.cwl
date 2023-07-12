#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Manta
doc: |-
  Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
  It is optimized for analysis of germline variation in small sets of individuals and somatic
  variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs,
  medium-sized indels and large insertions within a single efficient workflow. The method is
  designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is
  analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses
  can be completed within 2 hours. Manta combines paired and split-read evidence during SV
  discovery and scoring to improve accuracy, but does not require split-reads or successful
  breakpoint assemblies to report a variant in cases where there is strong evidence otherwise.

  It provides scoring models for germline variants in small sets of diploid samples and somatic
  variants in matched tumor/normal sample pairs. There is experimental support for analysis of
  unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and
  reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description
  of capabilities and limitations.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/manta:1.5.0

inputs:
- id: config
  label: config
  doc: |-
    provide a configuration file to override defaults in global config file (/opt/conda/share/manta-1.2.1-0/bin/configManta.py.ini)
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --config
    position: 1
    shellQuote: false
- id: bam
  label: bam
  doc: |-
    FILE Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [optional] (no default)
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: --bam
    position: 1
    shellQuote: false
- id: runDir
  label: runDir
  doc: |-
    Run script and run output will be written to this directory [required] (default: MantaWorkflow)
  type:
  - string
  - 'null'
  default: generated
  inputBinding:
    prefix: --runDir
    position: 1
    shellQuote: false
- id: reference
  label: reference
  doc: samtools-indexed reference fasta file [required]
  type: File
  secondaryFiles:
  - .fai
  inputBinding:
    prefix: --referenceFasta
    position: 1
    shellQuote: false
- id: tumorBam
  label: tumorBam
  doc: |-
    Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional=null]
  type:
  - File
  - 'null'
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: --tumorBam
    position: 1
    shellQuote: false
- id: exome
  label: exome
  doc: 'Set options for WES input: turn off depth filters'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --exome
    position: 1
    shellQuote: false
- id: rna
  label: rna
  doc: Set options for RNA-Seq input. Must specify exactly one bam input file
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --rna
    position: 1
    shellQuote: false
- id: unstrandedRNA
  label: unstrandedRNA
  doc: 'Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --unstrandedRNA
    position: 1
    shellQuote: false
- id: outputContig
  label: outputContig
  doc: Output assembled contig sequences in VCF file
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --outputContig
    position: 1
    shellQuote: false
- id: callRegions
  label: callRegions
  doc: |-
    Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
  inputBinding:
    prefix: --callRegions
    position: 1
    shellQuote: false
- id: mode
  label: mode
  doc: (-m) select run mode (local|sge)
  type: string
  default: local
  inputBinding:
    prefix: --mode
    position: 3
    shellQuote: false
- id: quiet
  label: quiet
  doc: |-
    Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --quiet
    position: 3
    shellQuote: false
- id: queue
  label: queue
  doc: (-q) specify scheduler queue name
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --queue
    position: 3
    shellQuote: false
- id: memgb
  label: memgb
  doc: |-
    (-g) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local  mode, 'unlimited' for sge mode)
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --memGb
    position: 3
    shellQuote: false
- id: maxTaskRuntime
  label: maxTaskRuntime
  doc: |-
    (format: hh:mm:ss) Specify scheduler max runtime per task, argument is provided to the 'h_rt' resource limit if using SGE (no default)
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --maxTaskRuntime
    position: 3
    shellQuote: false

outputs:
- id: python
  label: python
  type: File
  outputBinding:
    glob: $((inputs.runDir + "/runWorkflow.py"))
    outputEval: $((inputs.runDir.basename + "/runWorkflow.py"))
    loadContents: false
- id: pickle
  label: pickle
  type: File
  outputBinding:
    glob: $((inputs.runDir + "/runWorkflow.py.config.pickle"))
    outputEval: $((inputs.runDir.basename + "/runWorkflow.py.config.pickle"))
    loadContents: false
- id: candidateSV
  label: candidateSV
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.runDir + "/results/variants/candidateSV.vcf.gz"))
    outputEval: $((inputs.runDir.basename + "/results/variants/candidateSV.vcf.gz"))
    loadContents: false
- id: candidateSmallIndels
  label: candidateSmallIndels
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.runDir + "/results/variants/candidateSmallIndels.vcf.gz"))
    outputEval: $((inputs.runDir.basename + "/results/variants/candidateSmallIndels.vcf.gz"))
    loadContents: false
- id: diploidSV
  label: diploidSV
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.runDir + "/results/variants/diploidSV.vcf.gz"))
    outputEval: $((inputs.runDir.basename + "/results/variants/diploidSV.vcf.gz"))
    loadContents: false
- id: alignmentStatsSummary
  label: alignmentStatsSummary
  type: File
  outputBinding:
    glob: $((inputs.runDir + "/results/stats/alignmentStatsSummary.txt"))
    outputEval: $((inputs.runDir.basename + "/results/stats/alignmentStatsSummary.txt"))
    loadContents: false
- id: svCandidateGenerationStats
  label: svCandidateGenerationStats
  type: File
  outputBinding:
    glob: $((inputs.runDir + "/results/stats/svCandidateGenerationStats.tsv"))
    outputEval: $((inputs.runDir.basename + "/results/stats/svCandidateGenerationStats.tsv"))
    loadContents: false
- id: svLocusGraphStats
  label: svLocusGraphStats
  type: File
  outputBinding:
    glob: $((inputs.runDir + "/results/stats/svLocusGraphStats.tsv"))
    outputEval: $((inputs.runDir.basename + "/results/stats/svLocusGraphStats.tsv"))
    loadContents: false
- id: somaticSVs
  label: somaticSVs
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.runDir + "/results/variants/somaticSV.vcf.gz"))
    outputEval: $((inputs.runDir.basename + "/results/variants/somaticSV.vcf.gz"))
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments:
- position: 0
  valueFrom: configManta.py
  shellQuote: false
- position: 2
  valueFrom: $(";{runDir}/runWorkflow.py".replace(/\{runDir\}/g, inputs.runDir))
  shellQuote: false
- prefix: -j
  position: 3
  valueFrom: |-
    $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
  shellQuote: false
id: manta
