version development

task Gatk4SplitReads {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    String? outputFilename
    File bam
    File bam_bai
    File? intervals
    Array[String]? javaOptions
    Int? compression_level
    Boolean? addOutputSamProgramRecord
    Boolean? addOutputVcfCommandLine
    File? arguments_file
    String? cloudIndexPrefetchBuffer
    String? cloudPrefetchBuffer
    String? createOutputBamIndex
    String? createOutputBamMd5
    String? createOutputVariantIndex
    String? createOutputVariantMd5
    String? disableBamIndexCaching
    String? disableReadFilter
    Boolean? disableSequenceDictionaryValidation
    String? excludeIntervals
    File? gatkConfigFile
    Int? gcsRetries
    String? gcsProjectForRequesterPays
    Int? intervalExclusionPadding
    String? imr
    Int? ip
    String? isr
    Boolean? le
    Boolean? quiet
    String? readFilter
    String? readIndex
    String? readValidationStringency
    File? reference
    File? reference_fai
    File? reference_amb
    File? reference_ann
    File? reference_bwt
    File? reference_pac
    File? reference_sa
    File? reference_dict
    Float? secondsBetweenProgressUpdates
    String? sequenceDictionary
    Boolean? sitesOnlyVcfOutput
    String? splitLibraryName
    String? rg
    String? splitSample
    String? tmpDir
    Boolean? jdkDeflater
    Boolean? jdkInflater
    String? verbosity
    Boolean? disableToolDefaultReadFilters
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
    Boolean? invertSoftClipRatioFilter
    Float? softClippedLeadingTrailingRatio
    Float? softClippedRatioThreshold
  }
  command <<<
    set -e
    cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
    gatk SplitReads \
      --java-options '-Xmx~{((select_first([runtime_memory, 4, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
      --output '~{select_first([outputFilename, "."])}' \
      ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
      ~{if (defined(addOutputSamProgramRecord) && select_first([addOutputSamProgramRecord])) then "-add-output-sam-program-record" else ""} \
      ~{if (defined(addOutputVcfCommandLine) && select_first([addOutputVcfCommandLine])) then "-add-output-vcf-command-line" else ""} \
      ~{if defined(arguments_file) then ("--arguments_file:File '" + arguments_file + "'") else ""} \
      ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer '" + cloudIndexPrefetchBuffer + "'") else ""} \
      ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer '" + cloudPrefetchBuffer + "'") else ""} \
      ~{if defined(createOutputBamIndex) then ("--create-output-bam-index '" + createOutputBamIndex + "'") else ""} \
      ~{if defined(createOutputBamMd5) then ("--create-output-bam-md5 '" + createOutputBamMd5 + "'") else ""} \
      ~{if defined(createOutputVariantIndex) then ("--create-output-variant-index '" + createOutputVariantIndex + "'") else ""} \
      ~{if defined(createOutputVariantMd5) then ("--create-output-variant-md5 '" + createOutputVariantMd5 + "'") else ""} \
      ~{if defined(disableBamIndexCaching) then ("--disable-bam-index-caching '" + disableBamIndexCaching + "'") else ""} \
      ~{if defined(disableReadFilter) then ("--disable-read-filter '" + disableReadFilter + "'") else ""} \
      ~{if (defined(disableSequenceDictionaryValidation) && select_first([disableSequenceDictionaryValidation])) then "-disable-sequence-dictionary-validation" else ""} \
      ~{if defined(excludeIntervals) then ("--exclude-intervals '" + excludeIntervals + "'") else ""} \
      ~{if defined(gatkConfigFile) then ("--gatk-config-file '" + gatkConfigFile + "'") else ""} \
      ~{if defined(gcsRetries) then ("-gcs-retries " + gcsRetries) else ''} \
      ~{if defined(gcsProjectForRequesterPays) then ("--gcs-project-for-requester-pays '" + gcsProjectForRequesterPays + "'") else ""} \
      ~{if defined(intervalExclusionPadding) then ("--interval-exclusion-padding " + intervalExclusionPadding) else ''} \
      ~{if defined(imr) then ("-imr:IntervalMergingRule '" + imr + "'") else ""} \
      ~{if defined(ip) then ("-ip " + ip) else ''} \
      ~{if defined(isr) then ("-isr:IntervalSetRule '" + isr + "'") else ""} \
      ~{if (defined(le) && select_first([le])) then "--lenient" else ""} \
      ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
      ~{if defined(readFilter) then ("--read-filter '" + readFilter + "'") else ""} \
      ~{if defined(readIndex) then ("-read-index '" + readIndex + "'") else ""} \
      ~{if defined(readValidationStringency) then ("--read-validation-stringency '" + readValidationStringency + "'") else ""} \
      ~{if defined(reference) then ("--reference '" + reference + "'") else ""} \
      ~{if defined(secondsBetweenProgressUpdates) then ("-seconds-between-progress-updates " + secondsBetweenProgressUpdates) else ''} \
      ~{if defined(sequenceDictionary) then ("-sequence-dictionary '" + sequenceDictionary + "'") else ""} \
      ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output:Boolean" else ""} \
      ~{if defined(splitLibraryName) then ("--split-library-name '" + splitLibraryName + "'") else ""} \
      ~{if defined(rg) then ("--split-read-group '" + rg + "'") else ""} \
      ~{if defined(splitSample) then ("--split-sample '" + splitSample + "'") else ""} \
      ~{if defined(tmpDir) then ("--tmp-dir:GATKPathSpecifier '" + tmpDir + "'") else ""} \
      ~{if (defined(jdkDeflater) && select_first([jdkDeflater])) then "-jdk-deflater" else ""} \
      ~{if (defined(jdkInflater) && select_first([jdkInflater])) then "-jdk-inflater" else ""} \
      ~{if defined(verbosity) then ("-verbosity:LogLevel '" + verbosity + "'") else ""} \
      ~{if (defined(disableToolDefaultReadFilters) && select_first([disableToolDefaultReadFilters])) then "-disable-tool-default-read-filters" else ""} \
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
      ~{if defined(platformFilterName) then ("--platform-filter-name:String '" + platformFilterName + "'") else ""} \
      ~{if defined(blackListedLanes) then ("--black-listed-lanes:String '" + blackListedLanes + "'") else ""} \
      ~{if defined(readGroupBlackList) then ("--read-group-black-list:StringThe '" + readGroupBlackList + "'") else ""} \
      ~{if defined(keepReadGroup) then ("--keep-read-group:String '" + keepReadGroup + "'") else ""} \
      ~{if defined(maxReadLength) then ("--max-read-length " + maxReadLength) else ''} \
      ~{if defined(minReadLength) then ("--min-read-length " + minReadLength) else ''} \
      ~{if defined(readName) then ("--read-name:String '" + readName + "'") else ""} \
      ~{if (defined(keepReverseStrandOnly) && select_first([keepReverseStrandOnly])) then "--keep-reverse-strand-only" else ""} \
      ~{if defined(sample) then ("-sample:String '" + sample + "'") else ""} \
      ~{if (defined(invertSoftClipRatioFilter) && select_first([invertSoftClipRatioFilter])) then "--invert-soft-clip-ratio-filter" else ""} \
      ~{if defined(softClippedLeadingTrailingRatio) then ("--soft-clipped-leading-trailing-ratio " + softClippedLeadingTrailingRatio) else ''} \
      ~{if defined(softClippedRatioThreshold) then ("--soft-clipped-ratio-threshold " + softClippedRatioThreshold) else ''} \
      --input '~{bam}'
    if [ -f $(echo '~{basename(bam)}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{basename(bam)}' | sed 's/\.[^.]*$//').bai $(echo '~{basename(bam)}' ).bai; fi
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "broadinstitute/gatk:4.1.3.0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4, 4])}G"
    preemptible: 2
  }
  output {
    File out = basename(bam)
    File out_bai = basename(bam) + ".bai"
  }
}