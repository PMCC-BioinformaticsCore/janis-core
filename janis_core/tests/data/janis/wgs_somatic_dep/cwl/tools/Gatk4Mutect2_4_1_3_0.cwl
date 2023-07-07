#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: GatkMutect2
doc: |
  USAGE: Mutect2 [arguments]
  Call somatic SNVs and indels via local assembly of haplotypes
  Version:4.1.2.0

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
- id: tumorBams
  label: tumorBams
  doc: |-
    (--input) BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required. 
  type:
    type: array
    inputBinding:
      prefix: -I
    items: File
  inputBinding: {}
- id: normalBams
  label: normalBams
  doc: |-
    (--input) Extra BAM/SAM/CRAM file containing reads This argument must be specified at least once. Required. 
  type:
  - type: array
    inputBinding:
      prefix: -I
    items: File
  - 'null'
  inputBinding: {}
- id: normalSample
  label: normalSample
  doc: (--normal-sample, if) May be URL-encoded as output by GetSampleName with
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --normal-sample
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.vcf.gz
  inputBinding:
    prefix: -O
    position: 20
- id: reference
  label: reference
  doc: (-R) Reference sequence file Required.
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
    prefix: --reference
- id: outputBamName
  label: outputBamName
  doc: File to which assembled haplotypes should be written
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -bamout
- id: activityProfileOut
  label: activityProfileOut
  doc: 'Default value: null.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --activity-profile-out
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
- id: afOfAllelesNotInResource
  label: afOfAllelesNotInResource
  doc: |-
    (-default-af)  Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf fora derivation of the default value.  Default value: -1.0. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --af-of-alleles-not-in-resource
- id: alleles
  label: alleles
  doc: |-
    The set of alleles for which to force genotyping regardless of evidence Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --alleles
- id: annotation
  label: annotation
  doc: |-
    (-A) One or more specific annotations to add to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AlleleFraction, AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth, AS_ReadPosRankSumTest, AS_RMSMappingQuality, AS_StrandOddsRatio, BaseQuality, BaseQualityRankSumTest, ChromosomeCounts, ClippingRankSumTest, CountNs, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, ExcessHet, FisherStrand, FragmentLength, GenotypeSummaries, InbreedingCoeff, LikelihoodRankSumTest, MappingQuality, MappingQualityRankSumTest, MappingQualityZero, OrientationBiasReadCounts, OriginalAlignment, PossibleDeNovo, QualByDepth, ReadPosition, ReadPosRankSumTest, ReferenceBases, RMSMappingQuality, SampleList, StrandBiasBySample, StrandOddsRatio, TandemRepeat, UniqueAltReadCount}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --annotation
- id: annotationGroup
  label: annotationGroup
  doc: |-
    (-G) One or more groups of annotations to apply to variant calls This argument may be specified 0 or more times. Default value: null. Possible Values: {AS_StandardAnnotation, ReducibleAnnotation, StandardAnnotation, StandardHCAnnotation, StandardMutectAnnotation}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --annotation-group
- id: annotationsToExclude
  label: annotationsToExclude
  doc: |-
    (-AX)  One or more specific annotations to exclude from variant calls  This argument may be specified 0 or more times. Default value: null. Possible Values: {BaseQuality, Coverage, DepthPerAlleleBySample, DepthPerSampleHC, FragmentLength, MappingQuality, OrientationBiasReadCounts, ReadPosition, StrandBiasBySample, TandemRepeat}
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --annotations-to-exclude
- id: arguments_file
  label: arguments_file
  doc: |-
    read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. 
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --arguments_file
- id: assemblyRegionOut
  label: assemblyRegionOut
  doc: 'Output the assembly region to this IGV formatted file Default value: null.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --assembly-region-out
- id: baseQualityScoreThreshold
  label: baseQualityScoreThreshold
  doc: |2-
     Base qualities below this threshold will be reduced to the minimum (6)  Default value: 18.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --base-quality-score-threshold
- id: callableDepth
  label: callableDepth
  doc: |-
    Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. Default value: 10. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --callable-depth
- id: cloudIndexPrefetchBuffer
  label: cloudIndexPrefetchBuffer
  doc: |-
    (-CIPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --cloud-index-prefetch-buffer
- id: cloudPrefetchBuffer
  label: cloudPrefetchBuffer
  doc: |-
    (-CPB)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --cloud-prefetch-buffer
- id: createOutputBamIndex
  label: createOutputBamIndex
  doc: |-
    (-OBI)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --create-output-bam-index
- id: createOutputBamMd5
  label: createOutputBamMd5
  doc: |-
    (-OBM)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --create-output-bam-md5
- id: createOutputVariantIndex
  label: createOutputVariantIndex
  doc: |-
    (-OVI)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --create-output-variant-index
- id: createOutputVariantMd5
  label: createOutputVariantMd5
  doc: |-
    (-OVM)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --create-output-variant-md5
- id: disableBamIndexCaching
  label: disableBamIndexCaching
  doc: |-
    (-DBIC)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --disable-bam-index-caching
- id: disableReadFilter
  label: disableReadFilter
  doc: |-
    (-DF)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {GoodCigarReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotSecondaryAlignmentReadFilter, PassesVendorQualityCheckReadFilter, ReadLengthReadFilter, WellformedReadFilter}
  type:
  - boolean
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
- id: downsamplingStride
  label: downsamplingStride
  doc: |-
    (-stride)  Downsample a pool of reads starting within a range of one or more bases.  Default value: 1. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --downsampling-stride
- id: excludeIntervals
  label: excludeIntervals
  doc: '(-XLOne) This argument may be specified 0 or more times. Default value: null. '
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --exclude-intervals
- id: f1r2MaxDepth
  label: f1r2MaxDepth
  doc: 'sites with depth higher than this value will be grouped Default value: 200.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --f1r2-max-depth
- id: f1r2MedianMq
  label: f1r2MedianMq
  doc: 'skip sites with median mapping quality below this value Default value: 50.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --f1r2-median-mq
- id: f1r2MinBq
  label: f1r2MinBq
  doc: 'exclude bases below this quality from pileup Default value: 20.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --f1r2-min-bq
- id: f1r2TarGz_outputFilename
  label: f1r2TarGz_outputFilename
  doc: |-
    If specified, collect F1R2 counts and output files into this tar.gz file Default value: null. 
  type:
  - string
  - 'null'
  default: generated.tar.gz
  inputBinding:
    prefix: --f1r2-tar-gz
- id: founderId
  label: founderId
  doc: |-
    (--founder-id)  Samples representing the population founders This argument may be specified 0 or more times. Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -founder-id
- id: gatkConfigFile
  label: gatkConfigFile
  doc: 'A configuration file to use with the GATK. Default value: null.'
  type:
  - string
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
     Project to bill when accessing requester pays buckets. If unset, these buckets cannot be accessed.  Default value: . 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --gcs-project-for-requester-pays
- id: genotypeGermlineSites
  label: genotypeGermlineSites
  doc: |2-
     (EXPERIMENTAL) Call all apparent germline site even though they will ultimately be filtered.  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --genotype-germline-sites
- id: genotypePonSites
  label: genotypePonSites
  doc: |-
    Call sites in the PoN even though they will ultimately be filtered. Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --genotype-pon-sites
- id: germlineResource
  label: germlineResource
  doc: |2-
     Population vcf of germline sequencing containing allele fractions.  Default value: null. 
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
  inputBinding:
    prefix: --germline-resource
- id: graph
  label: graph
  doc: |-
    (--graph-output) Write debug assembly graph information to this file Default value: null.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -graph
- id: help
  label: help
  doc: |-
    (--help) display the help message Default value: false. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -h
- id: ignoreItrArtifacts
  label: ignoreItrArtifacts
  doc: ' inverted tandem repeats.  Default value: false. Possible values: {true, false} '
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --ignore-itr-artifactsTurn
- id: initialTumorLod
  label: initialTumorLod
  doc: |-
    (-init-lod)  Log 10 odds threshold to consider pileup active.  Default value: 2.0. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --initial-tumor-lod
- id: intervalExclusionPadding
  label: intervalExclusionPadding
  doc: |-
    (-ixp)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. 
  type:
  - string
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
    prefix: --interval-merging-rule
- id: ip
  label: ip
  doc: '(--interval-padding) Default value: 0.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -ipAmount
- id: isr
  label: isr
  doc: |-
    (--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --interval-set-rule
- id: intervals
  label: intervals
  doc: |-
    (-L) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. 
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --intervals
- id: le
  label: le
  doc: |-
    (--lenient) Lenient processing of VCF files Default value: false. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -LE
- id: maxPopulationAf
  label: maxPopulationAf
  doc: |-
    (-max-af)  Maximum population allele frequency in tumor-only mode.  Default value: 0.01. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --max-population-af
- id: maxReadsPerAlignmentStart
  label: maxReadsPerAlignmentStart
  doc: |2-
     Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.  Default value: 50. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-reads-per-alignment-start
- id: minBaseQualityScore
  label: minBaseQualityScore
  doc: |-
    (-mbq:Byte)  Minimum base quality required to consider a base for calling  Default value: 10. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --min-base-quality-score
- id: mitochondriaMode
  label: mitochondriaMode
  doc: |-
    Mitochondria mode sets emission and initial LODs to 0. Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --mitochondria-mode
- id: nativePairHmmThreads
  label: nativePairHmmThreads
  doc: ' How many threads should a native pairHMM implementation use  Default value:
    4. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --native-pair-hmm-threads
    valueFrom: |-
      $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
- id: nativePairHmmUseDoublePrecision
  label: nativePairHmmUseDoublePrecision
  doc: |2-
     use double precision in the native pairHmm. This is slower but matches the java implementation better  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --native-pair-hmm-use-double-precision
- id: normalLod
  label: normalLod
  doc: |-
    Log 10 odds threshold for calling normal variant non-germline. Default value: 2.2.
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --normal-lod
- id: encode
  label: encode
  doc: 'This argument may be specified 0 or more times. Default value: null.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -encode
- id: panelOfNormals
  label: panelOfNormals
  doc: |-
    (--panel-of-normals)  VCF file of sites observed in normal.  Default value: null. 
  type:
  - File
  - 'null'
  secondaryFiles:
  - .tbi
  inputBinding:
    prefix: --panel-of-normals
- id: pcrIndelQual
  label: pcrIndelQual
  doc: 'Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --pcr-indel-qual
- id: pcrSnvQual
  label: pcrSnvQual
  doc: 'Phred-scaled PCR SNV qual for overlapping fragments Default value: 40.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --pcr-snv-qual
- id: pedigree
  label: pedigree
  doc: |-
    (-ped) Pedigree file for determining the population founders. Default value: null.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --pedigree
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
    (-RF) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}
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
    (-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SILENT. Possible values: {STRICT, LENIENT, SILENT} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-validation-stringency
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
    prefix: --sites-only-vcf-output
- id: tmpDir
  label: tmpDir
  doc: 'Temp directory to use. Default value: null.'
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --tmp-dir
- id: tumorLodToEmit
  label: tumorLodToEmit
  doc: '(-emit-lod)  Log 10 odds threshold to emit variant to VCF.  Default value:
    3.0. '
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --tumor-lod-to-emit
- id: tumor
  label: tumor
  doc: |-
    (--tumor-sample) BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode argument.  Default value: null. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -tumor
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
    prefix: -verbosity
- id: version
  label: version
  doc: |-
    display the version number for this tool Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --version
- id: activeProbabilityThreshold
  label: activeProbabilityThreshold
  doc: |2-
     Minimum probability for a locus to be considered active.  Default value: 0.002. 
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --active-probability-threshold
- id: adaptivePruningInitialErrorRate
  label: adaptivePruningInitialErrorRate
  doc: ' Initial base error rate estimate for adaptive pruning  Default value: 0.001. '
  type:
  - double
  - 'null'
  inputBinding:
    prefix: --adaptive-pruning-initial-error-rate
- id: allowNonUniqueKmersInRef
  label: allowNonUniqueKmersInRef
  doc: |2-
     Allow graphs that have non-unique kmers in the reference  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --allow-non-unique-kmers-in-ref
- id: assemblyRegionPadding
  label: assemblyRegionPadding
  doc: |2-
     Number of additional bases of context to include around each assembly region  Default value: 100. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --assembly-region-padding
- id: bamWriterType
  label: bamWriterType
  doc: |-
    Which haplotypes should be written to the BAM Default value: CALLED_HAPLOTYPES. Possible values: {ALL_POSSIBLE_HAPLOTYPES, CALLED_HAPLOTYPES} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --bam-writer-type
- id: debugAssembly
  label: debugAssembly
  doc: |-
    (-debug)  Print out verbose debug information about each assembly region  Default value: false. Possible values: {true, false} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --debug-assembly
- id: disableAdaptivePruning
  label: disableAdaptivePruning
  doc: |2-
     Disable the adaptive algorithm for pruning paths in the graph  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --disable-adaptive-pruning
- id: disableToolDefaultAnnotations
  label: disableToolDefaultAnnotations
  doc: |-
    (--disable-tool-default-annotations)  Disable all tool default annotations  Default value: false. Possible values: {true, false}
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -disable-tool-default-annotations
- id: disableToolDefaultReadFilters
  label: disableToolDefaultReadFilters
  doc: |-
    (--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -disable-tool-default-read-filters
- id: dontIncreaseKmerSizesForCycles
  label: dontIncreaseKmerSizesForCycles
  doc: |2-
     Disable iterating over kmer sizes when graph cycles are detected  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --dont-increase-kmer-sizes-for-cycles
- id: dontTrimActiveRegions
  label: dontTrimActiveRegions
  doc: |2-
     If specified, we will not trim down the active region from the full region (active + extension) to just the active interval for genotyping  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --dont-trim-active-regions
- id: dontUseSoftClippedBases
  label: dontUseSoftClippedBases
  doc: |2-
     Do not analyze soft clipped bases in the reads  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --dont-use-soft-clipped-bases
- id: erc
  label: erc
  doc: |-
    (--emit-ref-confidence)  (BETA feature) Mode for emitting reference confidence scores  Default value: NONE. Possible values: {NONE, BP_RESOLUTION, GVCF} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -ERC
- id: enableAllAnnotations
  label: enableAllAnnotations
  doc: |2-
     Use all possible annotations (not for the faint of heart)  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --enable-all-annotations
- id: forceActive
  label: forceActive
  doc: |-
    If provided, all regions will be marked as active Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --force-active
- id: genotypeFilteredAlleles
  label: genotypeFilteredAlleles
  doc: |2-
     Whether to force genotype even filtered alleles  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --genotype-filtered-alleles
- id: gvcfLodBand
  label: gvcfLodBand
  doc: |-
    (-LODB) Exclusive upper bounds for reference confidence LOD bands (must be specified in increasing order)  This argument may be specified 0 or more times. Default value: [-2.5, -2.0, -1.5,
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --gvcf-lod-band
- id: kmerSize
  label: kmerSize
  doc: |-
    Kmer size to use in the read threading assembler This argument may be specified 0 or more times. Default value: [10, 25]. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --kmer-size
- id: maxAssemblyRegionSize
  label: maxAssemblyRegionSize
  doc: ' Maximum size of an assembly region  Default value: 300. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-assembly-region-size
- id: mnpDist
  label: mnpDist
  doc: |-
    (--max-mnp-distance)  Two or more phased substitutions separated by this distance or less are merged into MNPs.  Default value: 1. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -mnp-dist
- id: maxNumHaplotypesInPopulation
  label: maxNumHaplotypesInPopulation
  doc: |2-
     Maximum number of haplotypes to consider for your population  Default value: 128. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-num-haplotypes-in-population
- id: maxProbPropagationDistance
  label: maxProbPropagationDistance
  doc: |2-
     Upper limit on how many bases away probability mass can be moved around when calculating the boundaries between active and inactive assembly regions  Default value: 50. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-prob-propagation-distance
- id: maxSuspiciousReadsPerAlignmentStart
  label: maxSuspiciousReadsPerAlignmentStart
  doc: |2-
     Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.  Default value: 0. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-suspicious-reads-per-alignment-start
- id: maxUnprunedVariants
  label: maxUnprunedVariants
  doc: |2-
     Maximum number of variants in graph the adaptive pruner will allow  Default value: 100. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-unpruned-variants
- id: minAssemblyRegionSize
  label: minAssemblyRegionSize
  doc: ' Minimum size of an assembly region  Default value: 50. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-assembly-region-size
- id: minDanglingBranchLength
  label: minDanglingBranchLength
  doc: ' Minimum length of a dangling branch to attempt recovery  Default value: 4. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-dangling-branch-length
- id: minPruning
  label: minPruning
  doc: 'Minimum support to not prune paths in the graph Default value: 2.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --min-pruning
- id: minimumAlleleFraction
  label: minimumAlleleFraction
  doc: |-
    (-min-AF)  Lower bound of variant allele fractions to consider when calculating variant LOD  Default value: 0.0. 
  type:
  - float
  - 'null'
  inputBinding:
    prefix: --minimum-allele-fraction
- id: numPruningSamples
  label: numPruningSamples
  doc: 'Default value: 1.'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --num-pruning-samples
- id: pairHmmGapContinuationPenalty
  label: pairHmmGapContinuationPenalty
  doc: ' Flat gap continuation penalty for use in the Pair HMM  Default value: 10. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --pair-hmm-gap-continuation-penalty
- id: pairhmm
  label: pairhmm
  doc: |-
    (--pair-hmm-implementation)  The PairHMM implementation to use for genotype likelihood calculations  Default value: FASTEST_AVAILABLE. Possible values: {EXACT, ORIGINAL, LOGLESS_CACHING, AVX_LOGLESS_CACHING, AVX_LOGLESS_CACHING_OMP, EXPERIMENTAL_FPGA_LOGLESS_CACHING, FASTEST_AVAILABLE} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -pairHMM
- id: pcrIndelModel
  label: pcrIndelModel
  doc: |2-
     The PCR indel model to use  Default value: CONSERVATIVE. Possible values: {NONE, HOSTILE, AGGRESSIVE, CONSERVATIVE} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --pcr-indel-model
- id: phredScaledGlobalReadMismappingRate
  label: phredScaledGlobalReadMismappingRate
  doc: ' The global assumed mismapping rate for reads  Default value: 45. '
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --phred-scaled-global-read-mismapping-rate
- id: pruningLodThreshold
  label: pruningLodThreshold
  doc: 'Default value: 2.302585092994046. '
  type:
  - float
  - 'null'
  inputBinding:
    prefix: --pruning-lod-thresholdLn
- id: recoverAllDanglingBranches
  label: recoverAllDanglingBranches
  doc: |2-
     Recover all dangling branches  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: --recover-all-dangling-branches
- id: showhidden
  label: showhidden
  doc: |-
    (--showHidden)  display hidden arguments  Default value: false. Possible values: {true, false} 
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -showHidden
- id: smithWaterman
  label: smithWaterman
  doc: |2-
     Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice  Default value: JAVA. Possible values: {FASTEST_AVAILABLE, AVX_ENABLED, JAVA} 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --smith-waterman
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
    One or more genomic intervals to keep This argument must be specified at least once. Required. 
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --keep-intervals
- id: library
  label: library
  doc: |-
    (--library) Name of the library to keep This argument must be specified at least once. Required.
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
  doc: ' Minimum mapping quality to keep (inclusive)  Default value: 20. '
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
    prefix: --platform-filter-name
- id: blackListedLanes
  label: blackListedLanes
  doc: |-
    Platform unit (PU) to filter out This argument must be specified at least once. Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --black-listed-lanes
- id: readGroupBlackList
  label: readGroupBlackList
  doc: 'This argument must be specified at least once. Required. '
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --read-group-black-listThe
- id: keepReadGroup
  label: keepReadGroup
  doc: The name of the read group to keep Required.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: --keep-read-group
- id: maxReadLength
  label: maxReadLength
  doc: |-
    Keep only reads with length at most equal to the specified value Default value: 2147483647. 
  type:
  - int
  - 'null'
  inputBinding:
    prefix: --max-read-length
- id: minReadLength
  label: minReadLength
  doc: |-
    Keep only reads with length at least equal to the specified value Default value: 30.
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
    prefix: --read-name
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
    prefix: -sample

outputs:
- id: out
  label: out
  doc: To determine type
  type: File
  secondaryFiles:
  - .tbi
  outputBinding:
    glob: generated.vcf.gz
    loadContents: false
- id: stats
  label: stats
  doc: To determine type
  type: File
  outputBinding:
    glob: $((inputs.outputFilename + ".stats"))
    outputEval: $((inputs.outputFilename.basename + ".stats"))
    loadContents: false
- id: f1f2r_out
  label: f1f2r_out
  doc: To determine type
  type: File
  outputBinding:
    glob: generated.tar.gz
    loadContents: false
- id: bam
  label: bam
  doc: File to which assembled haplotypes should be written
  type:
  - File
  - 'null'
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
    glob: '$(inputs.outputBamName ? inputs.outputBamName : "generated")'
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand:
- gatk
- Mutect2
arguments:
- prefix: --java-options
  position: -1
  valueFrom: |-
    $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 16, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))
id: Gatk4Mutect2
