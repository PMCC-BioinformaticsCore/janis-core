version development

import "Gatk4MergeSamFiles_4_1_3_0.wdl" as G
import "Gatk4MarkDuplicates_4_1_3_0.wdl" as G2

workflow mergeAndMarkBams {
  input {
    Array[File] bams
    Array[File] bams_bai
    Boolean? createIndex = true
    Int? maxRecordsInRam = 5000000
    String? sampleName
    Boolean? mergeSamFiles_useThreading = true
    String? mergeSamFiles_validationStringency = "SILENT"
  }
  call G.Gatk4MergeSamFiles as mergeSamFiles {
    input:
      bams=bams,
      bams_bai=bams_bai,
      sampleName=sampleName,
      useThreading=select_first([mergeSamFiles_useThreading, true]),
      createIndex=select_first([createIndex, true]),
      maxRecordsInRam=select_first([maxRecordsInRam, 5000000]),
      validationStringency=select_first([mergeSamFiles_validationStringency, "SILENT"])
  }
  call G2.Gatk4MarkDuplicates as markDuplicates {
    input:
      bam=[mergeSamFiles.out],
      createIndex=select_first([createIndex, true]),
      maxRecordsInRam=select_first([maxRecordsInRam, 5000000])
  }
  output {
    File out = markDuplicates.out
    File out_bai = markDuplicates.out_bai
  }
}