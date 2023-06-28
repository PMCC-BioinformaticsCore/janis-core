cwlVersion: v1.0
class: CommandLineTool
id: structural_variants
label: structural_variants

doc: |
  Convert a VCF file to a BEDPE file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input)

hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/bioconductor-structuralvariantannotation:1.6.0--r40_0'

baseCommand: [Rscript]

inputs:
  input:
    type: File
    inputBinding:
      position: 2
  script:
    type: File
    inputBinding:
      position: 1
    default:
        class: File
        basename: "gridss.R"
        contents: |-
          library(stringr)
          library("VariantAnnotation")
          library("StructuralVariantAnnotation")

          args <- commandArgs(TRUE)
          vcfFile <- args[1]
          bedFile <- paste(sub('\\.vcf.gz$', '', vcfFile), ".bed", sep = "")

          simpleEventType <- function(gr) {
            return(
              ifelse(seqnames(gr) != seqnames(partner(gr)), "Translocation",
                     ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "Insertion",
                            ifelse(strand(gr) == strand(partner(gr)), "Inversion",
                                   ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                          "DUP")))))
          }
          options(scipen=999)
          full_vcf <- readVcf(vcfFile)
          gr <- breakpointRanges(full_vcf)
          bedpe <- data.frame(
            chrom1=seqnames(gr),
            start1=start(gr) - 1,
            end1=end(gr),
            chrom2=seqnames(partner(gr)),
            start2=start(partner(gr)) - 1,
            end2=end(partner(gr)),
            type=simpleEventType(gr),
            svLen=gr$svLen,
            insLen=gr$insLen,
            name=names(gr),
            score=gr$QUAL,
            strand1=strand(gr),
            strand2=strand(partner(gr))
          )

          bedpe <- bedpe[str_detect(bedpe$name, "gridss.+o"),]
          write.table(bedpe, bedFile, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.bed"
