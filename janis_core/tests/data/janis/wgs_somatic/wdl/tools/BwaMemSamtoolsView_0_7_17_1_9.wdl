version development

task BwaMemSamtoolsView {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File reference
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa
    File reference_dict
    Array[File] reads
    Array[File]? mates
    String? outputFilename
    String sampleName
    String? platformTechnology
    Int? minimumSeedLength
    Int? batchSize
    Boolean? useSoftClippingForSupplementaryAlignments
    Int? bandwidth
    Int? offDiagonalXDropoff
    Float? reseedTrigger
    Int? occurenceDiscard
    Boolean? performSW
    Int? matchingScore
    Int? mismatchPenalty
    Int? openGapPenalty
    Int? gapExtensionPenalty
    Int? clippingPenalty
    Int? unpairedReadPenalty
    Boolean? assumeInterleavedFirstInput
    Int? outputAlignmentThreshold
    Boolean? outputAllElements
    Boolean? appendComments
    Boolean? hardClipping
    Boolean? markShorterSplits
    Int? verboseLevel
    String? skippedReadsOutputFilename
    File? referenceIndex
    File? intervals
    String? includeReadsInReadGroup
    File? includeReadsInFile
    Int? includeReadsWithQuality
    String? includeReadsInLibrary
    Int? includeReadsWithCIGAROps
    Array[Int]? includeReadsWithAllFLAGs
    Array[Int]? includeReadsWithoutFLAGs
    Array[Int]? excludeReadsWithAllFLAGs
    Boolean? useMultiRegionIterator
    String? readTagToStrip
    Boolean? collapseBackwardCIGAROps
    String? outputFmt
  }
  command <<<
    set -e
     \
      bwa \
      mem \
      ~{reference} \
      ~{if defined(minimumSeedLength) then ("-k " + minimumSeedLength) else ''} \
      ~{if defined(batchSize) then ("-K " + batchSize) else ''} \
      ~{if (defined(useSoftClippingForSupplementaryAlignments) && select_first([useSoftClippingForSupplementaryAlignments])) then "-Y" else ""} \
      ~{if defined(bandwidth) then ("-w " + bandwidth) else ''} \
      ~{if defined(offDiagonalXDropoff) then ("-d " + offDiagonalXDropoff) else ''} \
      ~{if defined(reseedTrigger) then ("-r " + reseedTrigger) else ''} \
      ~{if defined(occurenceDiscard) then ("-c " + occurenceDiscard) else ''} \
      ~{if (defined(performSW) && select_first([performSW])) then "-P" else ""} \
      ~{if defined(matchingScore) then ("-A " + matchingScore) else ''} \
      ~{if defined(mismatchPenalty) then ("-B " + mismatchPenalty) else ''} \
      ~{if defined(openGapPenalty) then ("-O " + openGapPenalty) else ''} \
      ~{if defined(gapExtensionPenalty) then ("-E " + gapExtensionPenalty) else ''} \
      ~{if defined(clippingPenalty) then ("-L " + clippingPenalty) else ''} \
      ~{if defined(unpairedReadPenalty) then ("-U " + unpairedReadPenalty) else ''} \
      ~{if (defined(assumeInterleavedFirstInput) && select_first([assumeInterleavedFirstInput])) then "-p" else ""} \
      ~{if defined(outputAlignmentThreshold) then ("-T " + outputAlignmentThreshold) else ''} \
      ~{if (defined(outputAllElements) && select_first([outputAllElements])) then "-a" else ""} \
      ~{if (defined(appendComments) && select_first([appendComments])) then "-C" else ""} \
      ~{if (defined(hardClipping) && select_first([hardClipping])) then "-H" else ""} \
      ~{if (defined(markShorterSplits) && select_first([markShorterSplits])) then "-M" else ""} \
      ~{if defined(verboseLevel) then ("-v " + verboseLevel) else ''} \
      -R '@RG\tID:~{sampleName}\tSM:~{sampleName}\tLB:~{sampleName}\tPL:~{select_first([platformTechnology, "ILLUMINA"])}' \
      -t ~{select_first([runtime_cpu, 16, 1])} \
      ~{if length(reads) > 0 then sep(" ", reads) else ""} \
      ~{if (defined(mates) && length(select_first([mates])) > 0) then sep(" ", select_first([mates])) else ""} \
      | \
      samtools \
      view \
      -o ~{select_first([outputFilename, "~{sampleName}.bam"])} \
      ~{if defined(skippedReadsOutputFilename) then ("-U " + skippedReadsOutputFilename) else ''} \
      ~{if defined(referenceIndex) then ("-t " + referenceIndex) else ''} \
      ~{if defined(intervals) then ("-L " + intervals) else ''} \
      ~{if defined(includeReadsInReadGroup) then ("-r " + includeReadsInReadGroup) else ''} \
      ~{if defined(includeReadsInFile) then ("-R " + includeReadsInFile) else ''} \
      ~{if defined(includeReadsWithQuality) then ("-q " + includeReadsWithQuality) else ''} \
      ~{if defined(includeReadsInLibrary) then ("-l " + includeReadsInLibrary) else ''} \
      ~{if defined(includeReadsWithCIGAROps) then ("-m " + includeReadsWithCIGAROps) else ''} \
      ~{if (defined(includeReadsWithAllFLAGs) && length(select_first([includeReadsWithAllFLAGs])) > 0) then "-f " + sep(" ", select_first([includeReadsWithAllFLAGs])) else ""} \
      ~{if (defined(includeReadsWithoutFLAGs) && length(select_first([includeReadsWithoutFLAGs])) > 0) then "-F " + sep(" ", select_first([includeReadsWithoutFLAGs])) else ""} \
      ~{if (defined(excludeReadsWithAllFLAGs) && length(select_first([excludeReadsWithAllFLAGs])) > 0) then "-G " + sep(" ", select_first([excludeReadsWithAllFLAGs])) else ""} \
      ~{if (defined(useMultiRegionIterator) && select_first([useMultiRegionIterator])) then "-M" else ""} \
      ~{if defined(readTagToStrip) then ("-x " + readTagToStrip) else ''} \
      ~{if (defined(collapseBackwardCIGAROps) && select_first([collapseBackwardCIGAROps])) then "-B" else ""} \
      ~{if defined(outputFmt) then ("--output-fmt " + outputFmt) else ''} \
      -T ~{reference} \
      --threads ~{select_first([runtime_cpu, 16, 1])} \
      -h \
      -b
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 16, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "michaelfranklin/bwasamtools:0.7.17-1.9"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 16, 4])}G"
    preemptible: 2
  }
  output {
    File out = select_first([outputFilename, "~{sampleName}.bam"])
  }
}