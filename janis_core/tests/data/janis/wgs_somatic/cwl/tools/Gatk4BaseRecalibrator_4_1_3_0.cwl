#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: Base Recalibrator'
doc: |-
  First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
  The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

  This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
  that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
  We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
  Since there is a large amount of data one can then calculate an empirical probability of error given the 
  particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
  table (of the several covariate values, num observations, num mismatches, empirical quality score).

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.3.0

inputs:
- id: javaOptions
  label: javaOptions
  type:
  - type: array
    items: string
  - 'null'
- id: compression_level
  label: compression_level
  doc: |-
    Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
  type:
  - int
  - 'null'
- id: tmpDir
  label: tmpDir
  doc: Temp directory to use.
  type: string
  default: /tmp/
  inputBinding:
    prefix: --tmp-dir
- id: bam
  label: bam
  doc: BAM/SAM/CRAM file containing reads
  type: File
  secondaryFiles:
  - |-
    ${

            function resolveSecondary(base, secPattern) {
              if (secPattern[0] == "^") {
                var spl = base.split(".");
                var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
              }
              return base + secPattern
            }

            return [
                    {
                        location: resolveSecondary(self.location, "^.bai"),
                        basename: resolveSecondary(self.basename, ".bai"),
                        class: "File",
                    }
            ];

    }
  inputBinding:
    prefix: -I
    position: 6
- id: knownSites
  label: knownSites
  doc: |-
    **One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.** This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference, so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
  type:
    type: array
    inputBinding:
      prefix: --known-sites
    items: File
  inputBinding:
    position: 28
- id: reference
  label: reference
  doc: Reference sequence file
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: -R
    position: 5
- id: outputFilename
  label: outputFilename
  doc: |-
    **The output recalibration table filename to create.** After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate. Use '/dev/stdout' to print to standard out.
  type:
  - string
  - 'null'
  default: generated.table
  inputBinding:
    prefix: -O
    position: 8
    valueFrom: $(inputs.bam.basename).table
- id: intervals
  label: intervals
  doc: -L (BASE) One or more genomic intervals over which to operate
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --intervals
- id: intervalStrings
  label: intervalStrings
  doc: -L (BASE) One or more genomic intervals over which to operate
  type:
  - type: array
    inputBinding:
      prefix: --intervals
    items: string
  - 'null'
  inputBinding: {}

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: $(inputs.bam.basename).table
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- BaseRecalibrator
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 16, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4BaseRecalibrator
