version development

task Gatk4Mutect2 {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[String]? javaOptions
    Int? compression_level
    Array[File] tumorBams
    Array[File] tumorBams_bai
    Array[File]? normalBams
    Array[File]? normalBams_bai
    String? normalSample
    String? outputFilename
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    String? outputBamName
    String? activityProfileOut
    Boolean? addOutputSamProgramRecord
    Boolean? addOutputVcfCommandLine
    String? afOfAllelesNotInResource
    String? alleles
    String? annotation
    String? annotationGroup
    String? annotationsToExclude
    File? arguments_file
    String? assemblyRegionOut
    Int? baseQualityScoreThreshold
    Int? callableDepth
    Int? cloudIndexPrefetchBuffer
    Int? cloudPrefetchBuffer
    Boolean? createOutputBamIndex
    Boolean? createOutputBamMd5
    Boolean? createOutputVariantIndex
    Boolean? createOutputVariantMd5
    Boolean? disableBamIndexCaching
    Boolean? disableReadFilter
    Boolean? disableSequenceDictionaryValidation
    Int? downsamplingStride
    Boolean? excludeIntervals
    Int? f1r2MaxDepth
    Int? f1r2MedianMq
    Int? f1r2MinBq
    String? f1r2TarGz_outputFilename
    String? founderId
    String? gatkConfigFile
    Int? gcsRetries
    String? gcsProjectForRequesterPays
    Boolean? genotypeGermlineSites
    Boolean? genotypePonSites
    File? germlineResource
    File? germlineResource_tbi
    String? graph
    Boolean? help
    String? ignoreItrArtifacts
    String? initialTumorLod
    String? intervalExclusionPadding
    String? imr
    String? ip
    String? isr
    File? intervals
    Boolean? le
    String? maxPopulationAf
    Int? maxReadsPerAlignmentStart
    String? minBaseQualityScore
    Boolean? mitochondriaMode
    Int? nativePairHmmThreads
    Boolean? nativePairHmmUseDoublePrecision
    Float? normalLod
    String? encode
    File? panelOfNormals
    File? panelOfNormals_tbi
    Int? pcrIndelQual
    Int? pcrSnvQual
    String? pedigree
    Boolean? quiet
    String? readFilter
    String? readIndex
    String? readValidationStringency
    Float? secondsBetweenProgressUpdates
    String? sequenceDictionary
    Boolean? sitesOnlyVcfOutput
    String? tmpDir
    String? tumorLodToEmit
    String? tumor
    Boolean? jdkDeflater
    Boolean? jdkInflater
    String? verbosity
    Boolean? version
    Float? activeProbabilityThreshold
    Float? adaptivePruningInitialErrorRate
    Boolean? allowNonUniqueKmersInRef
    Int? assemblyRegionPadding
    String? bamWriterType
    String? debugAssembly
    Boolean? disableAdaptivePruning
    Boolean? disableToolDefaultAnnotations
    Boolean? disableToolDefaultReadFilters
    Boolean? dontIncreaseKmerSizesForCycles
    Boolean? dontTrimActiveRegions
    Boolean? dontUseSoftClippedBases
    String? erc
    Boolean? enableAllAnnotations
    Boolean? forceActive
    Boolean? genotypeFilteredAlleles
    String? gvcfLodBand
    Int? kmerSize
    Int? maxAssemblyRegionSize
    Int? mnpDist
    Int? maxNumHaplotypesInPopulation
    Int? maxProbPropagationDistance
    Int? maxSuspiciousReadsPerAlignmentStart
    Int? maxUnprunedVariants
    Int? minAssemblyRegionSize
    Int? minDanglingBranchLength
    Int? minPruning
    Float? minimumAlleleFraction
    Int? numPruningSamples
    Int? pairHmmGapContinuationPenalty
    String? pairhmm
    String? pcrIndelModel
    Int? phredScaledGlobalReadMismappingRate
    Float? pruningLodThreshold
    Boolean? recoverAllDanglingBranches
    Boolean? showhidden
    String? smithWaterman
    Int? ambigFilterBases
    Float? ambigFilterFrac
    Int? maxFragmentLength
    Int? minFragmentLength
    String? keepIntervals
    String? library
    Int? maximumMappingQuality
    Int? minimumMappingQuality
    Boolean? dontRequireSoftClipsBothEnds
    Int? filterTooShort
    String? platformFilterName
    String? blackListedLanes
    String? readGroupBlackList
    String? keepReadGroup
    Int? maxReadLength
    Int? minReadLength
    String? readName
    Boolean? keepReverseStrandOnly
    String? sample
  }
  command <<<
    set -e
    gatk Mutect2 \
      --java-options '-Xmx~{((select_first([runtime_memory, 16, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      ~{if length(tumorBams) > 0 then "-I '" + sep("' -I '", tumorBams) + "'" else ""} \
      ~{if (defined(normalBams) && length(select_first([normalBams])) > 0) then "-I '" + sep("' -I '", select_first([normalBams])) + "'" else ""} \
      ~{if defined(normalSample) then ("--normal-sample '" + normalSample + "'") else ""} \
      --reference '~{reference}' \
      ~{if defined(outputBamName) then ("-bamout '" + outputBamName + "'") else ""} \
      ~{if defined(activityProfileOut) then ("--activity-profile-out '" + activityProfileOut + "'") else ""} \
      ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "-add-output-sam-program-record" else ""} \
      ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "-add-output-vcf-command-line" else ""} \
      ~{if defined(afOfAllelesNotInResource) then ("--af-of-alleles-not-in-resource '" + afOfAllelesNotInResource + "'") else ""} \
      ~{if defined(alleles) then ("--alleles '" + alleles + "'") else ""} \
      ~{if defined(annotation) then ("--annotation '" + annotation + "'") else ""} \
      ~{if defined(annotationGroup) then ("--annotation-group '" + annotationGroup + "'") else ""} \
      ~{if defined(annotationsToExclude) then ("--annotations-to-exclude '" + annotationsToExclude + "'") else ""} \
      ~{if defined(arguments_file) then ("--arguments_file '" + arguments_file + "'") else ""} \
      ~{if defined(assemblyRegionOut) then ("--assembly-region-out '" + assemblyRegionOut + "'") else ""} \
      ~{if defined(baseQualityScoreThreshold) then ("--base-quality-score-threshold " + baseQualityScoreThreshold) else ''} \
      ~{if defined(callableDepth) then ("--callable-depth " + callableDepth) else ''} \
      ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer " + cloudIndexPrefetchBuffer) else ''} \
      ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer " + cloudPrefetchBuffer) else ''} \
      ~{if (defined(createOutputBamIndex) && select_first([createOutputBamIndex])) then "--create-output-bam-index" else ""} \
      ~{if (defined(createOutputBamMd5) && select_first([createOutputBamMd5])) then "--create-output-bam-md5" else ""} \
      ~{if (defined(createOutputVariantIndex) && select_first([createOutputVariantIndex])) then "--create-output-variant-index" else ""} \
      ~{if (defined(createOutputVariantMd5) && select_first([createOutputVariantMd5])) then "--create-output-variant-md5" else ""} \
      ~{if (defined(disableBamIndexCaching) && select_first([disableBamIndexCaching])) then "--disable-bam-index-caching" else ""} \
      ~{if (defined(disableReadFilter) && select_first([disableReadFilter])) then "--disable-read-filter" else ""} \
      ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "-disable-sequence-dictionary-validation" else ""} \
      ~{if defined(downsamplingStride) then ("--downsampling-stride " + downsamplingStride) else ''} \
      ~{if (defined(excludeIntervals) && select_first([excludeIntervals])) then "--exclude-intervals" else ""} \
      ~{if defined(f1r2MaxDepth) then ("--f1r2-max-depth " + f1r2MaxDepth) else ''} \
      ~{if defined(f1r2MedianMq) then ("--f1r2-median-mq " + f1r2MedianMq) else ''} \
      ~{if defined(f1r2MinBq) then ("--f1r2-min-bq " + f1r2MinBq) else ''} \
      --f1r2-tar-gz '~{select_first([f1r2TarGz_outputFilename, "generated.tar.gz"])}' \
      ~{if defined(founderId) then ("-founder-id '" + founderId + "'") else ""} \
      ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
      ~{if defined(gcsRetries) then ("-gcs-retries " + gcsRetries) else ''} \
      ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
      ~{if (defined(genotypeGermlineSites) && select_first([genotypeGermlineSites])) then "--genotype-germline-sites" else ""} \
      ~{if (defined(genotypePonSites) && select_first([genotypePonSites])) then "--genotype-pon-sites" else ""} \
      ~{if defined(germlineResource) then ("--germline-resource '" + germlineResource + "'") else ""} \
      ~{if defined(graph) then ("-graph '" + graph + "'") else ""} \
      ~{if (defined(help) && select_first([help])) then "-h" else ""} \
      ~{if defined(ignoreItrArtifacts) then ("--ignore-itr-artifactsTurn '" + ignoreItrArtifacts + "'") else ""} \
      ~{if defined(initialTumorLod) then ("--initial-tumor-lod '" + initialTumorLod + "'") else ""} \
      ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding '" + intervalExclusionPadding + "'") else ""} \
      ~{if defined(imr) then ("--interval-merging-rule '" + imr + "'") else ""} \
      ~{if defined(ip) then ("-ipAmount '" + ip + "'") else ""} \
      ~{if defined(isr) then ("--interval-set-rule '" + isr + "'") else ""} \
      ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
      ~{if (defined(le) && select_first([le])) then "-LE" else ""} \
      ~{if defined(maxPopulationAf) then ("--max-population-af '" + maxPopulationAf + "'") else ""} \
      ~{if defined(maxReadsPerAlignmentStart) then ("--max-reads-per-alignment-start " + maxReadsPerAlignmentStart) else ''} \
      ~{if defined(minBaseQualityScore) then ("--min-base-quality-score '" + minBaseQualityScore + "'") else ""} \
      ~{if (defined(mitochondriaMode) && select_first([mitochondriaMode])) then "--mitochondria-mode" else ""} \
      ~{if defined(select_first([nativePairHmmThreads, select_first([runtime_cpu, 1])])) then ("--native-pair-hmm-threads " + select_first([nativePairHmmThreads, select_first([runtime_cpu, 1])])) else ''} \
      ~{if (defined(nativePairHmmUseDoublePrecision) && select_first([nativePairHmmUseDoublePrecision])) then "--native-pair-hmm-use-double-precision" else ""} \
      ~{if defined(normalLod) then ("--normal-lod " + normalLod) else ''} \
      ~{if defined(encode) then ("-encode '" + encode + "'") else ""} \
      ~{if defined(panelOfNormals) then ("--panel-of-normals '" + panelOfNormals + "'") else ""} \
      ~{if defined(pcrIndelQual) then ("--pcr-indel-qual " + pcrIndelQual) else ''} \
      ~{if defined(pcrSnvQual) then ("--pcr-snv-qual " + pcrSnvQual) else ''} \
      ~{if defined(pedigree) then ("--pedigree '" + pedigree + "'") else ""} \
      ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
      ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
      ~{if defined(readIndex) then ("-read-index '" + readIndex + "'") else ""} \
      ~{if defined(readValidationStringency) then ("--read-validation-stringency '" + readValidationStringency + "'") else ""} \
      ~{if defined(secondsBetweenProgressUpdates) then ("-seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
      ~{if defined(sequenceDictionary) then ("-sequence-dictionary '" + sequenceDictionary + "'") else ""} \
      ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output" else ""} \
      ~{if defined(tmpDir) then ("--tmp-dir '" + tmpDir + "'") else ""} \
      ~{if defined(tumorLodToEmit) then ("--tumor-lod-to-emit '" + tumorLodToEmit + "'") else ""} \
      ~{if defined(tumor) then ("-tumor '" + tumor + "'") else ""} \
      ~{if (defined(jdkDeflater) && select_first([jdkDeflater])) then "-jdk-deflater" else ""} \
      ~{if (defined(jdkInflater) && select_first([jdkInflater])) then "-jdk-inflater" else ""} \
      ~{if defined(verbosity) then ("-verbosity '" + verbosity + "'") else ""} \
      ~{if (defined(version) && select_first([version])) then "--version" else ""} \
      ~{if defined(activeProbabilityThreshold) then ("--active-probability-threshold " + activeProbabilityThreshold) else ''} \
      ~{if defined(adaptivePruningInitialErrorRate) then ("--adaptive-pruning-initial-error-rate " + adaptivePruningInitialErrorRate) else ''} \
      ~{if (defined(allowNonUniqueKmersInRef) && select_first([allowNonUniqueKmersInRef])) then "--allow-non-unique-kmers-in-ref" else ""} \
      ~{if defined(assemblyRegionPadding) then ("--assembly-region-padding " + assemblyRegionPadding) else ''} \
      ~{if defined(bamWriterType) then ("--bam-writer-type '" + bamWriterType + "'") else ""} \
      ~{if defined(debugAssembly) then ("--debug-assembly '" + debugAssembly + "'") else ""} \
      ~{if (defined(disableAdaptivePruning) && select_first([disableAdaptivePruning])) then "--disable-adaptive-pruning" else ""} \
      ~{if (defined(disableToolDefaultAnnotations) && select_first([disableToolDefaultAnnotations])) then "-disable-tool-default-annotations" else ""} \
      ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "-disable-tool-default-read-filters" else ""} \
      ~{if (defined(dontIncreaseKmerSizesForCycles) && select_first([dontIncreaseKmerSizesForCycles])) then "--dont-increase-kmer-sizes-for-cycles" else ""} \
      ~{if (defined(dontTrimActiveRegions) && select_first([dontTrimActiveRegions])) then "--dont-trim-active-regions" else ""} \
      ~{if (defined(dontUseSoftClippedBases) && select_first([dontUseSoftClippedBases])) then "--dont-use-soft-clipped-bases" else ""} \
      ~{if defined(erc) then ("-ERC '" + erc + "'") else ""} \
      ~{if (defined(enableAllAnnotations) && select_first([enableAllAnnotations])) then "--enable-all-annotations" else ""} \
      ~{if (defined(forceActive) && select_first([forceActive])) then "--force-active" else ""} \
      ~{if (defined(genotypeFilteredAlleles) && select_first([genotypeFilteredAlleles])) then "--genotype-filtered-alleles" else ""} \
      ~{if defined(gvcfLodBand) then ("--gvcf-lod-band '" + gvcfLodBand + "'") else ""} \
      ~{if defined(kmerSize) then ("--kmer-size " + kmerSize) else ''} \
      ~{if defined(maxAssemblyRegionSize) then ("--max-assembly-region-size " + maxAssemblyRegionSize) else ''} \
      ~{if defined(mnpDist) then ("-mnp-dist " + mnpDist) else ''} \
      ~{if defined(maxNumHaplotypesInPopulation) then ("--max-num-haplotypes-in-population " + maxNumHaplotypesInPopulation) else ''} \
      ~{if defined(maxProbPropagationDistance) then ("--max-prob-propagation-distance " + maxProbPropagationDistance) else ''} \
      ~{if defined(maxSuspiciousReadsPerAlignmentStart) then ("--max-suspicious-reads-per-alignment-start " + maxSuspiciousReadsPerAlignmentStart) else ''} \
      ~{if defined(maxUnprunedVariants) then ("--max-unpruned-variants " + maxUnprunedVariants) else ''} \
      ~{if defined(minAssemblyRegionSize) then ("--min-assembly-region-size " + minAssemblyRegionSize) else ''} \
      ~{if defined(minDanglingBranchLength) then ("--min-dangling-branch-length " + minDanglingBranchLength) else ''} \
      ~{if defined(minPruning) then ("--min-pruning " + minPruning) else ''} \
      ~{if defined(minimumAlleleFraction) then ("--minimum-allele-fraction " + minimumAlleleFraction) else ''} \
      ~{if defined(numPruningSamples) then ("--num-pruning-samples " + numPruningSamples) else ''} \
      ~{if defined(pairHmmGapContinuationPenalty) then ("--pair-hmm-gap-continuation-penalty " + pairHmmGapContinuationPenalty) else ''} \
      ~{if defined(pairhmm) then ("-pairHMM '" + pairhmm + "'") else ""} \
      ~{if defined(pcrIndelModel) then ("--pcr-indel-model '" + pcrIndelModel + "'") else ""} \
      ~{if defined(phredScaledGlobalReadMismappingRate) then ("--phred-scaled-global-read-mismapping-rate " + phredScaledGlobalReadMismappingRate) else ''} \
      ~{if defined(pruningLodThreshold) then ("--pruning-lod-thresholdLn " + pruningLodThreshold) else ''} \
      ~{if (defined(recoverAllDanglingBranches) && select_first([recoverAllDanglingBranches])) then "--recover-all-dangling-branches" else ""} \
      ~{if (defined(showhidden) && select_first([showhidden])) then "-showHidden" else ""} \
      ~{if defined(smithWaterman) then ("--smith-waterman '" + smithWaterman + "'") else ""} \
      ~{if defined(ambigFilterBases) then ("--ambig-filter-bases " + ambigFilterBases) else ''} \
      ~{if defined(ambigFilterFrac) then ("--ambig-filter-frac " + ambigFilterFrac) else ''} \
      ~{if defined(maxFragmentLength) then ("--max-fragment-length " + maxFragmentLength) else ''} \
      ~{if defined(minFragmentLength) then ("--min-fragment-length " + minFragmentLength) else ''} \
      ~{if defined(keepIntervals) then ("--keep-intervals '" + keepIntervals + "'") else ""} \
      ~{if defined(library) then ("-library '" + library + "'") else ""} \
      ~{if defined(maximumMappingQuality) then ("--maximum-mapping-quality " + maximumMappingQuality) else ''} \
      ~{if defined(minimumMappingQuality) then ("--minimum-mapping-quality " + minimumMappingQuality) else ''} \
      ~{if (defined(dontRequireSoftClipsBothEnds) && select_first([dontRequireSoftClipsBothEnds])) then "--dont-require-soft-clips-both-ends" else ""} \
      ~{if defined(filterTooShort) then ("--filter-too-short " + filterTooShort) else ''} \
      ~{if defined(platformFilterName) then ("--platform-filter-name '" + platformFilterName + "'") else ""} \
      ~{if defined(blackListedLanes) then ("--black-listed-lanes '" + blackListedLanes + "'") else ""} \
      ~{if defined(readGroupBlackList) then ("--read-group-black-listThe '" + readGroupBlackList + "'") else ""} \
      ~{if defined(keepReadGroup) then ("--keep-read-group '" + keepReadGroup + "'") else ""} \
      ~{if defined(maxReadLength) then ("--max-read-length " + maxReadLength) else ''} \
      ~{if defined(minReadLength) then ("--min-read-length " + minReadLength) else ''} \
      ~{if defined(readName) then ("--read-name '" + readName + "'") else ""} \
      ~{if (defined(keepReverseStrandOnly) && select_first([keepReverseStrandOnly])) then "--keep-reverse-strand-only" else ""} \
      ~{if defined(sample) then ("-sample '" + sample + "'") else ""} \
      -O '~{select_first([outputFilename, "generated.vcf.gz"])}'
    if [ -f $(echo '~{outputBamName}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{outputBamName}' | sed 's/\.[^.]*$//').bai $(echo '~{outputBamName}' ).bai; fi
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 4, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.3.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 16, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "generated.vcf.gz"])
    File out_tbi = select_first([outputFilename, "generated.vcf.gz"]) + ".tbi"
    File stats = (select_first([outputFilename, "generated.vcf.gz"]) + ".stats")
    File f1f2r_out = select_first([f1r2TarGz_outputFilename, "generated.tar.gz"])
    File? bam = outputBamName
    File? bam_bai = if defined(outputBamName) then (outputBamName + ".bai") else None
  }
}