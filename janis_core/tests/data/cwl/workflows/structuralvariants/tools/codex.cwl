cwlVersion: v1.0
class: CommandLineTool
id: codex
label: codex

doc: |
  CODEX2

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)
      - $(inputs.mapping)
      - $(inputs.bed)

hints:
  DockerRequirement:
    dockerPull: 'migbro/codex2:3.8'

baseCommand: [Rscript]

inputs:
  input:
    type: 'File[]'
    inputBinding:
      position: 2
      itemSeparator: ','
  mapping:
    type: File
    inputBinding:
      position: 3
  bed:
    type: File
    inputBinding:
      position: 4
  chromosome:
    type: string
    inputBinding:
      position: 5
  script:
    type: File
    inputBinding:
      position: 1
    default:
        class: File
        basename: "CODEX2.R"
        contents: |-
          # load packages
          library(CODEX2)

          # parse arguments
          args <- commandArgs(TRUE)
          args2 <- args[2]

          bamsList <- as.vector(strsplit(args[1], ",")[[1]])
          bamsList <- sapply(strsplit(bamsList, "[/]"), "[[", 3)

          bams <- read.csv(args2, header = TRUE, sep = ",") # mapping file
          bedFile <- args[3]

          # mapping
          dirPath <- dirname(args2)
          bamsMap <- bams[bams$file %in% bamsList,]
          batch <- bamsMap$batch[[1]]
          batchname <- paste0("batch_", batch)

          # set global variables
          bamFile <- as.vector(bamsMap$file)
          bamdir <- file.path(dirPath, bamFile)
          sampname <- sapply(strsplit(bamFile, "[.]"), "[[", 1)

          bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, sampname = sampname, projectname = batchname)

          bamdir <- bambedObj$bamdir
          sampname <- bambedObj$sampname
          ref <- bambedObj$ref
          projectname <- bambedObj$projectname

          nsamples <- length(bamFile)

          ##########################################################
          # Getting GC content and mappability
          ##########################################################
          gc <- getgc(ref)
          mapp <- getmapp(ref)

          ##########################################################
          # Getting gene names, needed for targeted sequencing, here generating gene names in silico
          ##########################################################
          gene <- rep(NA, length(ref))
          for (chr in as.matrix(unique(seqnames(ref)))) {
            chr.index <- which(seqnames(ref) == chr)
            gene[chr.index] <- paste0(chr, "_gene_", ceiling(chr.index / 30))
          }
          values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))

          ##########################################################
          # Getting depth of coverage
          ##########################################################
          coverageObj <- getcoverage(bambedObj, mapqthres = 20)
          Y <- coverageObj$Y
          write.csv(Y, file = paste0(projectname, "_coverage.csv"), quote = FALSE)
          head(Y[, 1:nsamples])

          ##########################################################
          # Quality control
          ##########################################################
          qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, Inf), length_thresh = c(20, Inf), mapp_thresh = 0.9, gc_thresh = c(20, 80))
          Y_qc <- qcObj$Y_qc
          sampname_qc <- qcObj$sampname_qc
          ref_qc <- qcObj$ref_qc
          qcmat <- qcObj$qcmat
          gc_qc <- ref_qc$gc
          write.table(qcmat, file = paste0(projectname, "_qcmat", ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

          ##########################################################
          # Estimating library size factor for each sample
          ##########################################################
          Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x) { !any(x == 0) }), ]
          pseudo.sample <- apply(Y.nonzero, 1, function(x) { prod(x) ^ (1 / length(x)) })
          N <- apply(apply(Y.nonzero, 2, function(x) { x / pseudo.sample }), 2, median)
          plot(N, apply(Y, 2, sum), xlab = "Estimated library size factor", ylab = "Total sum of reads")

          ##########################################################
          # Genome-wide normalization using normalize_null
          ##########################################################
          # If there are negative control samples, use normalize_codex2_ns()
          # If there are negative control regions, use normalize_codex2_nr()
          normObj.null <- normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:nsamples, N = N)
          Yhat <- normObj.null$Yhat
          AIC <- normObj.null$AIC
          BIC <- normObj.null$BIC
          RSS <- normObj.null$RSS

          ##########################################################
          # Number of latent factors
          ##########################################################
          choiceofK(AIC, BIC, RSS, K = 1:nsamples , filename = "codex2_null_choiceofK.pdf")
          par(mfrow = c(1, 3))
          plot(1:nsamples, RSS, type = "b", xlab = "Number of latent variables", pch=20)
          plot(1:nsamples, AIC, type = "b", xlab = "Number of latent variables", pch=20)
          plot(1:nsamples, BIC, type = "b", xlab = "Number of latent variables", pch=20)
          par(mfrow = c(1,1))

          ##########################################################
          # CBS segmentation per chromosome: optimal for WGS and WES
          ##########################################################
          chr <- args[4] # FIXME: do this for each chromosome 
          optK <- which.max(BIC)
          chr.index <- which(seqnames(ref_qc) == chr)
          finalcall.CBS <- segmentCBS(Y_qc[chr.index,], Yhat, optK = which.max(BIC), K = 1:nsamples, sampname_qc = sampname_qc, ref_qc = ranges(ref_qc)[chr.index], chr = chr, lmax = 400, mode = "integer")

          # write results
          write.table(finalcall.CBS, file = paste0(projectname, "_", optK, ".CODEX_Chrom.txt"), sep = "\t", quote = F, row.names = F)

outputs:
  output:
    type: File
    outputBinding:
      glob: '*.CODEX_Chrom.txt'
