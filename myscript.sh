
successful builds:

LIMMA VOOM
    DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-limma:3.50.1--r41h5c21468_0' mulled-build --name-override 'janis-translate-limma-voom-3.34.9.9' --verbose build-and-test 'bioconductor-edger,r-statmod,r-scales,r-rjson,r-getopt,r-gplots,bioconductor-glimma'
    ->
    quay.io/biocontainers/janis-translate-limma-voom-3.34.9.9

HISAT2
    DEST_BASE_IMAGE='quay.io/biocontainers/hisat2:2.2.1--h87f3376_5' mulled-build --name-override 'janis-translate-hisat2-2.2.1' --verbose build-and-test 'samtools,seqtk'
    ->
    quay.io/biocontainers/janis-translate-hisat2-2.2.1

FEATURECOUNTS
    DEST_BASE_IMAGE='quay.io/biocontainers/coreutils:8.31--h14c3975_0' mulled-build --name-override 'janis-translate-featurecounts-2.0.1' --verbose build-and-test 'samtools,subread'
    ->
    quay.io/biocontainers/janis-translate-hisat2-2.2.1

# Run image [continuumio/miniconda3:latest] with command [[/bin/sh -c conda install  -c 'conda-forge' -c 'bioconda'  'bioconductor-edger' 'r-statmod' 'r-scales' 'r-rjson' 'r-getopt' 'r-gplots' 'bioconductor-glimma' --strict-channel-priority -p /usr/local --copy --yes --quiet]]

FGSEA
    DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-fgsea:1.24.0--r42hc247a5b_0' mulled-build --name-override 'janis-translate-fgsea-1.24.0' --verbose build-and-test 'r-optparse'
    ->
    quay.io/biocontainers/janis-translate-fgsea-1.8.0 

EGSEA
    DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-egsea:1.26.0--r42hdfd78af_0' mulled-build --name-override 'janis-translate-egsea-1.26.0' --verbose build-and-test 'r-optparse,r-rjson,r-statmod'
    ->
    quay.io/biocontainers/janis-translate-egsea-1.26.0

GOSEQ
    DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-goseq:1.44.0--r41hdfd78af_0' mulled-build --name-override 'janis-translate-goseq-1.44.0' --verbose build-and-test 'bioconductor-goseq,bioconductor-org.hs.eg.db,bioconductor-org.dm.eg.db,bioconductor-org.dr.eg.db,bioconductor-org.mm.eg.db,r-dplyr,r-ggplot2,r-optparse'
    ->
    quay.io/biocontainers/janis-translate-goseq-1.44.0



mulled-build --verbose build-and-test 'quay.io/biocontainers/r-statmod:1.4.29--r3.3.2_0,quay.io/biocontainers/r-scales:0.4.1--r3.3.2_1,quay.io/biocontainers/r-rjson:0.2.15--r3.3.2_0,quay.io/biocontainers/r-getopt:1.20.0--r3.3.2_0,quay.io/biocontainers/r-gplots:2.17.0--r3.3.2_0'

DEST_BASE_IMAGE='quay.io/biocontainers/r-statmod:1.4.29--r3.3.2_0' mulled-build --verbose build-and-test 'r-scales=0.4.1--r3.3.2_1,r-rjson=0.2.15--r3.3.2_0,r-getopt=1.20.0--r3.3.2_0,r-gplots=2.17.0--r3.3.2_0'


bioconductor-edger=3.20.7
r-statmod=1.4.30
r-scales=0.5.0
r-rjson=0.2.15
r-getopt=1.20.0
r-gplots=3.0.1
bioconductor-glimma=1.6.0

DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-limma:3.50.1--r41h5c21468_0' mulled-build --verbose build-and-test 'bioconductor-edger=3.20.7,r-statmod=1.4.30,r-scales=0.5.0,r-rjson=0.2.15,r-getopt=1.20.0,r-gplots=3.0.1,bioconductor-glimma=1.6.0'

DEST_BASE_IMAGE='quay.io/biocontainers/bioconductor-limma:3.50.1--r41h5c21468_0' mulled-build --verbose build-and-test 'bioconductor-edger,r-statmod,r-scales,r-rjson,r-getopt,r-gplots,bioconductor-glimma'


# quay.io/biocontainers/bioconductor-limma:3.50.1--r41h5c21468_0
# quay.io/biocontainers/bioconductor-edger:3.36.0--r41h619a076_1 
# quay.io/biocontainers/r-statmod:1.4.29--r3.3.2_0
# quay.io/biocontainers/r-scales:0.4.1--r3.3.2_1
# quay.io/biocontainers/r-rjson:0.2.15--r3.3.2_0
# quay.io/biocontainers/r-getopt:1.20.0--r3.3.2_0
# quay.io/biocontainers/r-gplots:2.17.0--r3.3.2_0
# quay.io/biocontainers/bioconductor-glimma:2.4.0--r41hdfd78af_0