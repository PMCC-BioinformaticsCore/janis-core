#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Strelka (Somatic)
doc: |-
  Usage: configureStrelkaSomaticWorkflow.py [options]
  Version: 2.9.10
  This script configures Strelka somatic small variant calling.
  You must specify an alignment file (BAM or CRAM) for each sample of a matched tumor-normal pair.
  Configuration will produce a workflow run script which can execute the workflow on a single node or through
  sge and resume any interrupted execution.

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/strelka:2.9.10

inputs:
- id: normalBam
  label: normalBam
  doc: Normal sample BAM or CRAM file. (no default)
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: --normalBam=
    position: 1
    separate: false
- id: tumorBam
  label: tumorBam
  doc: (--tumorBam)  Tumor sample BAM or CRAM file. [required] (no default)
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: --tumourBam=
    position: 1
    separate: false
- id: reference
  label: reference
  doc: ' samtools-indexed reference fasta file [required]'
  type: File
  secondaryFiles:
  - .fai
  inputBinding:
    prefix: --referenceFasta=
    position: 1
    separate: false
- id: rundir
  label: rundir
  doc: |-
    Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory. (default: StrelkaSomaticWorkflow)
  type:
  - string
  - 'null'
  default: generated
  inputBinding:
    prefix: --runDir=
    position: 1
    separate: false
- id: region
  label: region
  doc: |-
    Limit the analysis to one or more genome region(s) for debugging purposes. If this argument is provided multiple times the union of all specified regions will be analyzed. All regions must be non-overlapping to get a meaningful result. Examples: '--region chr20' (whole chromosome), '--region chr2:100-2000 --region chr3:2500-3000' (two regions)'. If this option is specified (one or more times) together with the 'callRegions' BED file,then all region arguments will be intersected with the callRegions BED track.
  type:
  - type: array
    inputBinding:
      prefix: --region
    items: string
  - 'null'
  inputBinding:
    position: 1
- id: config
  label: config
  doc: |-
    provide a configuration file to override defaults in global config file (/opt/strelka/bin/configureStrelkaSomaticWorkflow.py.ini)
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --config=
    position: 1
    separate: false
- id: outputcallableregions
  label: outputcallableregions
  doc: Output a bed file describing somatic callable regions of the genome
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --outputCallableRegions
    position: 1
    separate: true
- id: indelCandidates
  label: indelCandidates
  doc: |-
    Specify a VCF of candidate indel alleles. These alleles are always evaluated but only reported in the output when they are inferred to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored. This option may be specified more than once, multiple input VCFs will be merged. (default: None)
  type:
  - type: array
    inputBinding:
      prefix: --indelCandidates=
      separate: false
    items: File
  - 'null'
  inputBinding:
    position: 1
- id: forcedgt
  label: forcedgt
  doc: |-
    Specify a VCF of candidate alleles. These alleles are always evaluated and reported even if they are unlikely to exist in the sample. The VCF must be tabix indexed. All indel alleles must be left- shifted/normalized, any unnormalized allele will trigger a runtime error. This option may be specified more than once, multiple input VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF, excluded from block compression), but the specific SNV alleles are ignored. (default: None)
  type:
  - type: array
    inputBinding:
      prefix: --forcedGT=
      separate: false
    items: File
  - 'null'
  inputBinding:
    position: 1
- id: targeted
  label: targeted
  doc: |-
    Set options for other targeted input: note in particular that this flag turns off high-depth filters
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --targeted
    position: 1
    separate: true
- id: exome
  label: exome
  doc: |-
    Set options for exome: note in particular that this flag turns off high-depth filters
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --exome
    position: 1
    separate: true
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
    prefix: --callRegions=
    position: 1
    separate: false
- id: noisevcf
  label: noisevcf
  doc: Noise vcf file (submit argument multiple times for more than one file)
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
  inputBinding:
    prefix: --noiseVcf=
    position: 1
    separate: false
- id: scansizemb
  label: scansizemb
  doc: |-
    Maximum sequence region size (in megabases) scanned by each task during genome variant calling. (default: 12)
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --scanSizeMb=
    position: 1
    separate: false
- id: callmemmb
  label: callmemmb
  doc: |-
    Set variant calling task memory limit (in megabytes). It is not recommended to change the default in most cases, but this might be required for a sample of unusual depth.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --callMemMb=
    position: 1
    separate: false
- id: retaintempfiles
  label: retaintempfiles
  doc: Keep all temporary files (for workflow debugging)
  type: boolean
  default: false
  inputBinding:
    prefix: --retainTempFiles
    position: 1
    separate: true
- id: disableevs
  label: disableevs
  doc: Disable empirical variant scoring (EVS).
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --disableEVS
    position: 1
    separate: true
- id: reportevsfeatures
  label: reportevsfeatures
  doc: ' Report all empirical variant scoring features in VCF output.'
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --reportEVSFeatures
    position: 1
    separate: true
- id: snvscoringmodelfile
  label: snvscoringmodelfile
  doc: |2-
     Provide a custom empirical scoring model file for SNVs (default: /opt/strelka/share/config/somaticSNVScoringM odels.json)
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --snvScoringModelFile=
    position: 1
    separate: false
- id: indelscoringmodelfile
  label: indelscoringmodelfile
  doc: |2-
     Provide a custom empirical scoring model file for indels (default: /opt/strelka/share/config/somaticInde lScoringModels.json)
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --indelScoringModelFile=
    position: 1
    separate: false
- id: mode
  label: mode
  doc: (-m MODE)  select run mode (local|sge)
  type: string
  default: local
  inputBinding:
    prefix: --mode
    position: 3
    shellQuote: false
- id: queue
  label: queue
  doc: (-q QUEUE) specify scheduler queue name
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --queue
    position: 3
    shellQuote: false
- id: memGb
  label: memGb
  doc: |2-
     (-g MEMGB) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local mode, 'unlimited' for sge mode)
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --memGb
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

outputs:
- id: configPickle
  label: configPickle
  type: File
  outputBinding:
    glob: $((inputs.rundir + "/runWorkflow.py.config.pickle"))
    outputEval: $((inputs.rundir.basename + "/runWorkflow.py.config.pickle"))
    loadContents: false
- id: script
  label: script
  type: File
  outputBinding:
    glob: $((inputs.rundir + "/runWorkflow.py"))
    outputEval: $((inputs.rundir.basename + "/runWorkflow.py"))
    loadContents: false
- id: stats
  label: stats
  doc: |-
    A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
  type: File
  outputBinding:
    glob: $((inputs.rundir + "/results/stats/runStats.tsv"))
    outputEval: $((inputs.rundir.basename + "/results/stats/runStats.tsv"))
    loadContents: false
- id: indels
  label: indels
  doc: ''
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.rundir + "/results/variants/somatic.indels.vcf.gz"))
    outputEval: $((inputs.rundir.basename + "/results/variants/somatic.indels.vcf.gz"))
    loadContents: false
- id: snvs
  label: snvs
  doc: ''
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: $((inputs.rundir + "/results/variants/somatic.snvs.vcf.gz"))
    outputEval: $((inputs.rundir.basename + "/results/variants/somatic.snvs.vcf.gz"))
    loadContents: false
stdout: _stdout
stderr: _stderr
arguments:
- position: 0
  valueFrom: configureStrelkaSomaticWorkflow.py
- position: 2
  valueFrom: $(";{rundir}/runWorkflow.py".replace(/\{rundir\}/g, inputs.rundir))
  shellQuote: false
- prefix: --jobs
  position: 3
  valueFrom: |-
    $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
  shellQuote: false
id: strelka_somatic
