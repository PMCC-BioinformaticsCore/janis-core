#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: 'GATK4: SplitReads'
doc: |-
  USAGE: SplitReads [arguments]
  Outputs reads from a SAM/BAM/CRAM by read group, sample and library name
  Version:4.1.3.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.1.3.0

inputs:
- id: outputFilename
  label: outputFilename
  doc: "The directory to output SAM/BAM/CRAM files. Default value: '.' "
  type: string
  default: .
  inputBinding:
    prefix: --output
- id: bam
  label: bam
  doc: |-
    (-I:String) BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
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
    prefix: --input
    position: 1
- id: intervals
  label: intervals
  doc: |-
    (-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --intervals
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
- id: addOutputSamProgramRecord
  label: addOutputSamProgramRecord
  doc: |-
    (--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -add-output-sam-program-record
- id: addOutputVcfCommandLine
  label: addOutputVcfCommandLine
  doc: |-
    (--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -add-output-vcf-command-line
- id: arguments_file
  label: arguments_file
  doc: |-
    read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --arguments_file:File
- id: cloudIndexPrefetchBuffer
  label: cloudIndexPrefetchBuffer
  doc: |-
    (-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --cloud-index-prefetch-buffer
- id: cloudPrefetchBuffer
  label: cloudPrefetchBuffer
  doc: |-
    (-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --cloud-prefetch-buffer
- id: createOutputBamIndex
  label: createOutputBamIndex
  doc: |-
    (-OBI:Boolean)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --create-output-bam-index
- id: createOutputBamMd5
  label: createOutputBamMd5
  doc: |-
    (-OBM:Boolean)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --create-output-bam-md5
- id: createOutputVariantIndex
  label: createOutputVariantIndex
  doc: |-
    (-OVI:Boolean)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --create-output-variant-index
- id: createOutputVariantMd5
  label: createOutputVariantMd5
  doc: |-
    (-OVM:Boolean)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --create-output-variant-md5
- id: disableBamIndexCaching
  label: disableBamIndexCaching
  doc: |-
    (-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --disable-bam-index-caching
- id: disableReadFilter
  label: disableReadFilter
  doc: |-
    (-DF:String)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --disable-read-filter
- id: disableSequenceDictionaryValidation
  label: disableSequenceDictionaryValidation
  doc: |-
    (--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -disable-sequence-dictionary-validation
- id: excludeIntervals
  label: excludeIntervals
  doc: |-
    (-XL:StringOne) This argument may be specified 0 or more times. Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --exclude-intervals
- id: gatkConfigFile
  label: gatkConfigFile
  doc: 'A configuration file to use with the GATK. Default value: null.'
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --gatk-config-file
- id: gcsRetries
  label: gcsRetries
  doc: |-
    (--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -gcs-retries
- id: gcsProjectForRequesterPays
  label: gcsProjectForRequesterPays
  doc: |2-
     Project to bill when accessing requester pays  buckets. If unset, these buckets cannot be accessed.  Default value: . 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --gcs-project-for-requester-pays
- id: intervalExclusionPadding
  label: intervalExclusionPadding
  doc: |-
    (-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --interval-exclusion-padding
- id: imr
  label: imr
  doc: |-
    (--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -imr:IntervalMergingRule
- id: ip
  label: ip
  doc: '(--interval-padding) Default value: 0.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -ip
- id: isr
  label: isr
  doc: |-
    (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -isr:IntervalSetRule
- id: le
  label: le
  doc: |-
    (-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --lenient
- id: quiet
  label: quiet
  doc: |-
    Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --QUIET
- id: readFilter
  label: readFilter
  doc: |-
    (-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-filter
- id: readIndex
  label: readIndex
  doc: |-
    (--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -read-index
- id: readValidationStringency
  label: readValidationStringency
  doc: |-
    (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SITool returned: 0 LENT. Possible values: {STRICT, LENIENT, SILENT} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-validation-stringency
- id: reference
  label: reference
  doc: '(-R:String) Reference sequence Default value: null.'
  type:
  - File
  - 'null'
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: --reference
- id: secondsBetweenProgressUpdates
  label: secondsBetweenProgressUpdates
  doc: |-
    (--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0. 
  type:
  - double
  - 'null'
  inputBinding:
    prefix: -seconds-between-progress-updates
- id: sequenceDictionary
  label: sequenceDictionary
  doc: |-
    (--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -sequence-dictionary
- id: sitesOnlyVcfOutput
  label: sitesOnlyVcfOutput
  doc: |2-
     If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --sites-only-vcf-output:Boolean
- id: splitLibraryName
  label: splitLibraryName
  doc: |-
    (-LB)  Split file by library.  Default value: false. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --split-library-name
- id: rg
  label: rg
  doc: '(-RG:BooleanSplit) Default value: false. Possible values: {true, false}'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --split-read-group
- id: splitSample
  label: splitSample
  doc: |-
    (-SM:Boolean) Split file by sample. Default value: false. Possible values: {true, false}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --split-sample
- id: tmpDir
  label: tmpDir
  doc: 'Temp directory to use. Default value: null.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --tmp-dir:GATKPathSpecifier
- id: jdkDeflater
  label: jdkDeflater
  doc: |-
    (--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -jdk-deflater
- id: jdkInflater
  label: jdkInflater
  doc: |-
    (--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -jdk-inflater
- id: verbosity
  label: verbosity
  doc: |-
    (--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -verbosity:LogLevel
- id: disableToolDefaultReadFilters
  label: disableToolDefaultReadFilters
  doc: |-
    (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -disable-tool-default-read-filters
- id: ambigFilterBases
  label: ambigFilterBases
  doc: |-
    Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --ambig-filter-bases
- id: ambigFilterFrac
  label: ambigFilterFrac
  doc: |-
    Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --ambig-filter-frac
- id: maxFragmentLength
  label: maxFragmentLength
  doc: 'Default value: 1000000.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-fragment-length
- id: minFragmentLength
  label: minFragmentLength
  doc: 'Default value: 0.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-fragment-length
- id: keepIntervals
  label: keepIntervals
  doc: |-
    Valid only if "IntervalOverlapReadFilter" is specified: One or more genomic intervals to keep This argument must be specified at least once. Required. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --keep-intervals
- id: library
  label: library
  doc: |-
    (--library) Valid only if "LibraryReadFilter" is specified: Name of the library to keep This argument must be specified at least once. Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -library
- id: maximumMappingQuality
  label: maximumMappingQuality
  doc: ' Maximum mapping quality to keep (inclusive)  Default value: null. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --maximum-mapping-quality
- id: minimumMappingQuality
  label: minimumMappingQuality
  doc: ' Minimum mapping quality to keep (inclusive)  Default value: 10. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --minimum-mapping-quality
- id: dontRequireSoftClipsBothEnds
  label: dontRequireSoftClipsBothEnds
  doc: |2-
     Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --dont-require-soft-clips-both-ends
- id: filterTooShort
  label: filterTooShort
  doc: 'Minimum number of aligned bases Default value: 30.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --filter-too-short
- id: platformFilterName
  label: platformFilterName
  doc: This argument must be specified at least once. Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --platform-filter-name:String
- id: blackListedLanes
  label: blackListedLanes
  doc: |-
    Platform unit (PU) to filter out This argument must be specified at least once. Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --black-listed-lanes:String
- id: readGroupBlackList
  label: readGroupBlackList
  doc: 'This argument must be specified at least once. Required. '
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-group-black-list:StringThe
- id: keepReadGroup
  label: keepReadGroup
  doc: The name of the read group to keep Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --keep-read-group:String
- id: maxReadLength
  label: maxReadLength
  doc: Keep only reads with length at most equal to the specified value Required.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-read-length
- id: minReadLength
  label: minReadLength
  doc: |-
    Keep only reads with length at least equal to the specified value Default value: 1.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-read-length
- id: readName
  label: readName
  doc: Keep only reads with this read name Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-name:String
- id: keepReverseStrandOnly
  label: keepReverseStrandOnly
  doc: |2-
     Keep only reads on the reverse strand  Required. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --keep-reverse-strand-only
- id: sample
  label: sample
  doc: |-
    (--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -sample:String
- id: invertSoftClipRatioFilter
  label: invertSoftClipRatioFilter
  doc: |2-
     Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --invert-soft-clip-ratio-filter
- id: softClippedLeadingTrailingRatio
  label: softClippedLeadingTrailingRatio
  doc: |2-
     Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --soft-clipped-leading-trailing-ratio
- id: softClippedRatioThreshold
  label: softClippedRatioThreshold
  doc: |2-
     Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --soft-clipped-ratio-threshold

outputs:
- id: out
  label: out
  doc: Bam
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
                        path: resolveSecondary(self.path, "^.bai"),
                        basename: resolveSecondary(self.basename, ".bai"),
                        class: "File",
                    }
            ];

    }
  outputBinding:
    glob: $(inputs.bam.basename)
    outputEval: $(inputs.bam.basename.basename)
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- SplitReads
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4SplitReads
