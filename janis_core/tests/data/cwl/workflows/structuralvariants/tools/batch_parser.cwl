cwlVersion: v1.0
class: CommandLineTool
id: batch_parser
label: batch_parser

doc: |
  Group samples by batch

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)
      - $(inputs.samples)

hints:
  DockerRequirement:
    dockerPull: 'rocker/r-base:4.0.5'

baseCommand: [Rscript]

inputs:
  input:
    type: 'File[]'
    inputBinding:
      position: 2
      itemSeparator: ','
  samples:
    type: File
    inputBinding:
      position: 3
  script:
    type: File
    inputBinding:
      position: 1
    default:
        class: File
        basename: "batch_parser.R"
        contents: |-
          args <- commandArgs(TRUE)

          bamsList <- as.vector(strsplit(args[1], ",")[[1]])
          bamsList <- sapply(strsplit(bamsList, "[/]"), "[[", 3)
          samplesFile <- read.delim(args[2], header = T, sep = " ")

          bams <- data.frame(unlist(bamsList))
          colnames(bams)[1] <- "file"
          samples <- samplesFile

          indx <- sapply(samples$sample_id, grep, bams$file)
          temp_df <- cbind(samples[unlist(indx), , drop = F], bams[unlist(indx), , drop = F])
          temp_df$index_file <- paste(temp_df$file, ".bai", sep = "")
          output <- split(temp_df, temp_df$batch)

          sapply(names(output), function(x)
            write.table(output[[x]], file = paste0("bams_x_batch_", x, ".csv"), sep = ",", row.names = F))

outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: 'bams_x_batch_*.csv'
