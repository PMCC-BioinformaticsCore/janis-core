
import os
import operator
from datetime import date
from abc import ABC, abstractmethod

from janis_core import (
    ToolInput,
    Filename,
    File,
    Int,
    String,
    Boolean,
    ToolOutput,
    InputSelector,
    ToolMetadata,
    ToolArgument,
    Float,
    Array
)

from janis_core.tool.test_classes import (
    TTestPreprocessor,
    TTestExpectedOutput,
    TTestCase,
)

from ..types import TextFile, BamBai, Bam, FastaWithDict, Sam, Cram
from .bioinformaticstool import BioinformaticsTool
from janis_core.types import UnionType




### BASE ###

class SamTools_1_7:
    def container(self):
        return "biocontainers/samtools:v1.7.0_cv3"

    def version(self):
        return "1.7.0"

class SamTools_1_9:
    def container(self):
        return "quay.io/biocontainers/samtools:1.9--h8571acd_11"

    def version(self):
        return "1.9.0"
    

class SamToolsToolBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "Samtools"

    @classmethod
    @abstractmethod
    def samtools_command(cls):
        raise Exception(
            "Subclass must implement the samtools_command method: expects one of: ["
            "   view, sort, index, idxstats, flagstat, stats, bedcov, depth, "
            "   merge, faidx, fqidx, tview, split, quickcheck, dict, fixmate, "
            "   mpileup, flags, fastq/a, collate, refheader, cat, rmdup, "
            "   addreplacerg, calmd, targetcut, phase, depad, markdup"
            "]"
        )

    @classmethod
    def base_command(cls):
        return ["samtools", cls.samtools_command()]

    def inputs(self):
        return []

    def doc(self):
        return """
    Samtools is a set of utilities that manipulate alignments in the BAM format. It imports from 
    and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, 
    and allows to retrieve reads in any regions swiftly.

    Samtools is designed to work on a stream. It regards an input file `-' as the standard input (stdin) 
    and an output file `-' as the standard output (stdout). Several commands can thus be combined with 
    Unix pipes. Samtools always output warning and error messages to the standard error output (stderr).

    Samtools is also able to open a BAM (not SAM) file on a remote FTP or HTTP server if the BAM file 
    name starts with `ftp://' or `http://'. Samtools checks the current working directory for the index 
    file and will download the index upon absence. Samtools does not retrieve the entire alignment file 
    unless it is asked to do so.

    Documentation: http://www.htslib.org/doc/samtools.html#DESCRIPTION""".strip()

    @abstractmethod
    def container(self):
        raise Exception(
            "An error likely occurred when resolving the method order for docker for the samtools classes "
            "or you're trying to execute the docker method of the base class (ie, don't do that). "
            "The method order resolution must preference Gatkbase subclasses, "
            "and the subclass must contain a definition for docker."
        )

    def arguments(self):
        return []


### VIEW ###

class SamToolsViewBase(SamToolsToolBase, ABC):
    def tool(self):
        return "SamToolsView"

    @classmethod
    def samtools_command(cls):
        return "view"

    def inputs(self):
        return [
            *super(SamToolsViewBase, self).inputs(),
            *SamToolsViewBase.additional_inputs,
            ToolInput("sam", UnionType(Sam(), Bam(), Cram()), position=10),
            ToolInput(
                "reference",
                FastaWithDict(optional=True),
                position=6,
                prefix="-T",
                doc="A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed "
                "by samtools faidx. If an index is not present, one will be generated for you.",
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("sam", remove_file_extension=True),
                    extension=".bam",
                ),
                position=5,
                prefix="-o",
                doc="Output to FILE [stdout].",
            ),
            ToolInput(
                "regions",
                Array(String, optional=True),
                position=11,
                doc="Region specifications after the input filename to restrict output to only those alignments which "
                "overlap the specified region(s). Use of region specifications requires a coordinate-sorted and "
                "indexed input file (in BAM or CRAM format)",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Bam(), glob=InputSelector("outputFilename"))]

    def friendly_name(self):
        return "SamTools: View"

    def bind_metadata(self):
        self.metadata = ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Samtools",
            doi=None,
            citation=None,  # find citation
            keywords=["samtools", "view"],
            documentationUrl="http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS",
            documentation="""Ensure SAMTOOLS.SORT is inheriting from parent metadata
        
---------------------------------------------------------------------------------------------------
    
With no options or regions specified, prints all alignments in the specified input alignment file 
(in SAM, BAM, or CRAM format) to standard output in SAM format (with no header).

You may specify one or more space-separated region specifications after the input filename to 
restrict output to only those alignments which overlap the specified region(s). 
Use of region specifications requires a coordinate-sorted and indexed input file (in BAM or CRAM format).""".strip(),
        )
        return self.metadata

    def arguments(self):
        return [
            ToolArgument(
                "-S",
                position=2,
                doc="Ignored for compatibility with previous samtools versions. Previously this option was "
                "required if input was in SAM format, but now the correct format is automatically "
                "detected by examining the first few characters of input.",
            ),
            ToolArgument("-h", position=3, doc="Include the header in the output."),
            ToolArgument("-b", position=4, doc="Output in the BAM format."),
        ]

    additional_inputs = [
        ToolInput(
            "cramOutput",
            Boolean(optional=True),
            position=5,
            prefix="-C",
            doc="Output in the CRAM format (requires -T).",
        ),
        ToolInput(
            "compressedBam",
            Boolean(optional=True),
            position=5,
            prefix="-1",
            doc="Enable fast BAM compression (implies -b).",
        ),
        ToolInput(
            "uncompressedBam",
            Boolean(optional=True),
            position=5,
            prefix="-u",
            doc="Output uncompressed BAM. This option saves time spent on compression/decompression and is "
            "thus preferred when the output is piped to another samtools command.",
        ),
        ToolInput(
            "onlyOutputHeader",
            Boolean(optional=True),
            position=5,
            prefix="-H",
            doc="Output the header only.",
        ),
        ToolInput(
            "countAlignments",
            Boolean(optional=True),
            position=5,
            prefix="-c",
            doc="Instead of printing the alignments, only count them and print the total number. "
            "All filter options, such as -f, -F, and -q, are taken into account.",
        ),
        # ToolInput("", Boolean(), position=5, prefix="-?", doc="Output long help and exit immediately."),
        ToolInput(
            "writeAlignments",
            File(optional=True),
            position=5,
            prefix="-U",
            doc="Write alignments that are not selected by the various filter options to FILE. "
            "When this option is used, all alignments (or all alignments intersecting the regions specified) "
            "are written to either the output file or this file, but never both.",
        ),
        ToolInput(
            "inputTSV",
            File(optional=True),
            position=5,
            prefix="-t",
            doc="A tab-delimited FILE. Each line must contain the reference name in the first column and the "
            "length of the reference in the second column, with one line for each distinct reference. "
            "Any additional fields beyond the second column are ignored. This file also defines the order "
            "of the reference sequences in sorting. If you run: `samtools faidx <ref.fa>', the resulting "
            "index file <ref.fa>.fai can be used as this FILE.",
        ),
        ToolInput(
            "onlyOverlapping",
            File(optional=True),
            position=5,
            prefix="-L",
            doc="Only output alignments overlapping the input BED FILE [null].",
        ),
        ToolInput(
            "useMultiRegionIterator",
            Boolean(optional=True),
            position=5,
            prefix="-M",
            doc="Use the multi-region iterator on the union of the BED file and command-line region arguments. "
            "This avoids re-reading the same regions of files so can sometimes be much faster. "
            "Note this also removes duplicate sequences. Without this a sequence that overlaps multiple "
            "regions specified on the command line will be reported multiple times.",
        ),
        ToolInput(
            "outputAlignmentsInReadGroup",
            String(optional=True),
            position=5,
            prefix="-r",
            doc="Output alignments in read group STR [null]. Note that records with no RG tag will also be "
            "output when using this option. This behaviour may change in a future release.",
        ),
        ToolInput(
            "outputAlignmentsInFileReadGroups",
            File(optional=True),
            position=5,
            prefix="-R",
            doc="Output alignments in read groups listed in FILE [null]. Note that records with no RG tag "
            "will also be output when using this option. This behaviour may change in a future release.",
        ),
        ToolInput(
            "mapqThreshold",
            Int(optional=True),
            position=5,
            prefix="-q",
            doc="Skip alignments with MAPQ smaller than INT [0].",
        ),
        ToolInput(
            "outputAlignmentsInLibrary",
            String(optional=True),
            position=5,
            prefix="-l",
            doc="Only output alignments in library STR [null].",
        ),
        ToolInput(
            "outputAlignmentsMeetingCIGARThreshold",
            Int(optional=True),
            position=5,
            prefix="-m",
            doc="Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]",
        ),
        ToolInput(
            "outputAlignmentsWithBitsSet",
            String(optional=True),
            position=5,
            prefix="-f",
            doc="Only output alignments with all bits set in INT present in the FLAG field. "
            "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
            "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].",
        ),
        ToolInput(
            "doNotOutputAlignmentsWithBitsSet",
            String(optional=True),
            position=5,
            prefix="-F",
            doc="Do not output alignments with any bits set in INT present in the FLAG field. "
            "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
            "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].",
        ),
        ToolInput(
            "doNotOutputAlignmentsWithAllBitsSet",
            String(optional=True),
            position=5,
            prefix="-G",
            doc="Do not output alignments with all bits set in INT present in the FLAG field. "
            "This is the opposite of -f such that -f12 -G12 is the same as no filtering at all. "
            "INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or "
            "in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].",
        ),
        ToolInput(
            "readTagToExclude",
            String(optional=True),
            position=5,
            prefix="-x",
            doc="Read tag to exclude from output (repeatable) [null]",
        ),
        ToolInput(
            "collapseBackwardCIGAR",
            Boolean(optional=True),
            position=5,
            prefix="-B",
            doc="Collapse the backward CIGAR operation.",
        ),
        ToolInput(
            "subsamplingProportion",
            Float(optional=True),
            position=5,
            prefix="-s",
            doc="Output only a proportion of the input alignments. This subsampling acts in the same "
            "way on all of the alignment records in the same template or read pair, so it never "
            "keeps a read but not its mate. The integer and fractional parts of the -s INT.FRAC "
            "option are used separately: the part after the decimal point sets the fraction of "
            "templates/pairs to be kept, while the integer part is used as a seed that influences "
            "which subset of reads is kept.",
        ),
        ToolInput(
            "threads",
            Int(optional=True),
            position=5,
            prefix="-@",
            doc="Number of BAM compression threads to use in addition to main thread [0].",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "sam": f"{remote_dir}/NA12878-BRCA1.bwamem.stdout",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "threads": 16,
                },
                output=Bam.basic_test(
                    "out",
                    2740774,
                    f"{remote_dir}/NA12878-BRCA1.bam.flagstat",
                    "9a6af420f287df52a122ac723f41b535",
                ),
            ),
            TTestCase(
                name="minimal",
                input={
                    "sam": f"{remote_dir}/NA12878-BRCA1.bwamem.stdout",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "threads": 16,
                },
                output=self.minimal_test(),
            ),
        ]


class SamToolsView_1_7(SamTools_1_7, SamToolsViewBase):
    pass


class SamToolsView_1_9(SamTools_1_9, SamToolsViewBase):
    pass


SamToolsViewLatest = SamToolsView_1_9





### FLAGSTAT ###

class SamToolsFlagstatBase(SamToolsToolBase, ABC):
    def tool(self):
        return "SamToolsFlagstat"

    @classmethod
    def samtools_command(cls):
        return "flagstat"

    def inputs(self):
        return [
            ToolInput("bam", Bam(), position=10),
            ToolInput(
                "threads",
                Int(optional=True),
                position=5,
                prefix="-@",
                doc="Number of BAM compression threads to use in addition to main thread [0].",
            ),
            ToolInput("outputFilename", Filename, prefix=">", position=11),
        ]

    def outputs(self):
        return [ToolOutput("out", TextFile, selector=InputSelector("outputFilename"))]

    def friendly_name(self):
        return "SamTools: Flagstat"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Jiaan Yu"],
            dateCreated=date(2020, 2, 14),
            dateUpdated=date(2020, 2, 14),
            institution="Samtools",
            doi=None,
            citation=None,
            keywords=["samtools", "flagstat"],
            documentationUrl="http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS",
            documentation="""Does a full pass through the input file to calculate and print statistics to stdout.

Provides counts for each of 13 categories based primarily on bit flags in the FLAG field. Each category in the output is broken down into QC pass and QC fail. In the default output format, these are presented as "#PASS + #FAIL" followed by a description of the category.

The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:

122 + 28 in total (QC-passed reads + QC-failed reads)

Which would indicate that there are a total of 150 reads in the input file, 122 of which are marked as QC pass and 28 of which are marked as "not passing quality controls"

Following this, additional categories are given for reads which are:

secondary     0x100 bit set

supplementary     0x800 bit set

duplicates     0x400 bit set

mapped     0x4 bit not set

paired in sequencing     0x1 bit set

read1     both 0x1 and 0x40 bits set

read2     both 0x1 and 0x80 bits set

properly paired     both 0x1 and 0x2 bits set and 0x4 bit not set

with itself and mate mapped     0x1 bit set and neither 0x4 nor 0x8 bits set

singletons     both 0x1 and 0x8 bits set and bit 0x4 not set

And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:

with mate mapped to a different chr     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME

with mate mapped to a different chr (mapQ>=5)     0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5)""".strip(),
        )

    # def tests(self):
    #     return [
    #         TTestCase(
    #             name="basic",
    #             input={
    #                 "bam": os.path.join(
    #                     BioinformaticsTool.test_data_path(), "small.bam"
    #                 ),
    #             },
    #             output=[
    #                 TTestExpectedOutput(
    #                     tag="out",
    #                     preprocessor=TTestPreprocessor.FileMd5,
    #                     operator=operator.eq,
    #                     expected_value="dc58fe92a9bb0c897c85804758dfadbf",
    #                 ),
    #                 TTestExpectedOutput(
    #                     tag="out",
    #                     preprocessor=TTestPreprocessor.FileContent,
    #                     operator=operator.contains,
    #                     expected_value="19384 + 0 in total (QC-passed reads + QC-failed reads)",
    #                 ),
    #                 TTestExpectedOutput(
    #                     tag="out",
    #                     preprocessor=TTestPreprocessor.LineCount,
    #                     operator=operator.eq,
    #                     expected_value=13,
    #                 ),
    #             ],
    #         )
    #     ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                },
                output=TextFile.basic_test(
                    "out",
                    410,
                    "19486 + 0 in total (QC-passed reads + QC-failed reads)",
                    13,
                    "ddbcfe52e60b925d222fb8bc1517a7a0",
                ),
            )
        ]


class SamToolsFlagstat_1_7(SamTools_1_7, SamToolsFlagstatBase):
    pass


class SamToolsFlagstat_1_9(SamTools_1_9, SamToolsFlagstatBase):
    pass


SamToolsFlagstatLatest = SamToolsFlagstat_1_9




### MPILEUP ###

class SamToolsMpileupBase(SamToolsToolBase, ABC):
    def tool(self):
        return "SamToolsMpileup"

    @classmethod
    def samtools_command(cls):
        return "mpileup"

    def inputs(self):
        return [
            *self.additional_inputs,
            ToolInput("bam", BamBai(), position=10),
        ]

    def outputs(self):
        return [ToolOutput("out", TextFile, glob=InputSelector("outputFilename"))]

    def friendly_name(self):
        return "SamTools: Mpileup"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Jiaan Yu"],
            dateCreated=date(2020, 5, 19),
            dateUpdated=date(2020, 5, 19),
            institution="Samtools",
            doi=None,
            citation=None,
            keywords=["samtools", "mpileup"],
            documentationUrl="http://www.htslib.org/doc/samtools-mpileup.html",
            documentation="""Generate text pileup output for one or multiple BAM files. Each input file produces a separate group of pileup columns in the output.

Samtools mpileup can still produce VCF and BCF output (with -g or -u), but this feature is deprecated and will be removed in a future release. Please use bcftools mpileup for this instead. (Documentation on the deprecated options has been removed from this manual page, but older versions are available online at <http://www.htslib.org/doc/>.)

Note that there are two orthogonal ways to specify locations in the input file; via -r region and -l file. The former uses (and requires) an index to do random access while the latter streams through the file contents filtering out the specified regions, requiring no index. The two may be used in conjunction. For example a BED file containing locations of genes in chromosome 20 could be specified using -r 20 -l chr20.bed, meaning that the index is used to find chromosome 20 and then it is filtered for the regions listed in the bed file.""".strip(),
        )

    additional_inputs = [
        ToolInput(
            "illuminaEncoding",
            Boolean(optional=True),
            prefix="--illumina1.3+",
            doc="Assume the quality is in the Illumina 1.3+ encoding.",
        ),
        ToolInput(
            "countOrphans",
            Boolean(optional=True),
            prefix="--count-orphans",
            doc="do not discard anomalous read pairs",
        ),
        # Not sure this would load the
        # ToolInput("bamList", File(optional=True), prefix="--bam-list", doc="list of input BAM filenames, one per line")
        ToolInput(
            "noBAQ",
            Boolean(optional=True),
            prefix="--no-BAQ",
            doc="disable BAQ (per-Base Alignment Quality)",
        ),
        ToolInput(
            "adjustMQ",
            Int(optional=True),
            prefix="--adjust-MQ",
            doc="adjust mapping quality; recommended:50, disable:0 [0]",
        ),
        ToolInput(
            "maxDepth",
            Int(optional=True),
            prefix="--max-depth",
            doc="max per-file depth; avoids excessive memory usage [8000]",
        ),
        ToolInput(
            "redoBAQ",
            Boolean(optional=True),
            prefix="--redo-BAQ",
            doc="recalculate BAQ on the fly, ignore existing BQs",
        ),
        ToolInput(
            "fastaRef",
            File(optional=True),
            prefix="--fasta-ref",
            doc=" skip unlisted positions (chr pos) or regions (BED)",
        ),
        ToolInput(
            "excludeRG",
            File(optional=True),
            prefix="--exclude-RG",
            doc="exclude read groups listed in FILE",
        ),
        ToolInput(
            "positions",
            File(optional=True),
            prefix="--positions",
            doc="skip unlisted positions (chr pos) or regions (BED)",
        ),
        ToolInput(
            "minBQ",
            Int(optional=True),
            prefix="--min-BQ",
            doc="Minimum base quality for a base to be considered [13]",
        ),
        ToolInput(
            "minMQ",
            Int(optional=True),
            prefix="--min-MQ",
            doc="skip alignments with mapQ smaller than INT [0]",
        ),
        ToolInput(
            "region",
            String(optional=True),
            prefix="--region",
            doc="region in which pileup is generated",
        ),
        ToolInput(
            "ignoreRG",
            Boolean(optional=True),
            prefix="--ignore-RG",
            doc="ignore RG tags (one BAM = one sample)",
        ),
        ToolInput(
            "inclFlags",
            String(optional=True),
            prefix="--incl-flags",
            doc="required flags: skip reads with mask bits unset []",
        ),
        ToolInput(
            "exclFlags",
            String(optional=True),
            prefix="--excl-flags",
            doc="filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]",
        ),
        ToolInput(
            "outputFilename",
            Filename(extension=".txt"),
            prefix="--output",
            doc="write output to FILE [standard output]",
        ),
        ToolInput(
            "ignoreOverlaps",
            Boolean(optional=True),
            prefix="--ignore-overlaps",
            doc="disable read-pair overlap detection",
        ),
        ToolInput(
            "outputBP",
            Boolean(optional=True),
            prefix="--output-BP",
            doc="output base positions on reads",
        ),
        ToolInput(
            "outputMQ",
            Boolean(optional=True),
            prefix="--output-MQ",
            doc="output mapping quality",
        ),
        ToolInput(
            "outputQNAME",
            Boolean(optional=True),
            prefix="--output-QNAME",
            doc="output read names",
        ),
        ToolInput(
            "allPositions",
            Boolean(optional=True),
            prefix="-a",
            doc="output all positions (including zero depth)",
        ),
        ToolInput(
            "absolutelyAllPositions",
            Boolean(optional=True),
            doc="output absolutely all positions, including unused ref. sequences",
        ),
        ToolInput(
            "reference",
            FastaWithDict(optional=True),
            prefix="--reference",
            doc="Reference sequence FASTA FILE [null]",
        ),
    ]

    """
    def tests(self):
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": os.path.join(
                        BioinformaticsTool.test_data_path(), "small.bam"
                    ),
                },
                output=[
                    TTestExpectedOutput(
                        tag="out",
                        preprocessor=TTestPreprocessor.FileMd5,
                        operator=operator.eq,
                        expected_value="6b6f2401df9965b5250f4752dde03f2a",
                    ),
                    TTestExpectedOutput(
                        tag="out",
                        preprocessor=TTestPreprocessor.FileContent,
                        operator=operator.contains,
                        expected_value="17:43044045-43125733\t5\tN\t15\tCCCCCCCCCCCCCCC\tJDDAJDEDCDJD>gB\n",
                    ),
                    TTestExpectedOutput(
                        tag="out",
                        preprocessor=TTestPreprocessor.LineCount,
                        operator=operator.eq,
                        expected_value=81689,
                    ),
                ],
            )
        ]
        """

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "positions": f"{remote_dir}/NA12878-BRCA1.sorted.uncompressed.stdout",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "countOrphans": True,
                    "noBAQ": True,
                    "maxDepth": 10000,
                    "minBQ": 0,
                },
                output=TextFile.basic_test(
                    "out",
                    19900,
                    "chr17\t43044391\tG\t19\tA,A,,A.a,,A,,A..,,a\tDJCJ:FHDDBJBBJJIDDB",
                    187,
                    "53c3e03c20730ff45411087444379b1b",
                ),
            )
        ]



class SamToolsMpileup_1_7(SamTools_1_7, SamToolsMpileupBase):
    pass


class SamToolsMpileup_1_9(SamTools_1_9, SamToolsMpileupBase):
    pass


SamToolsMpileupLatest = SamToolsMpileup_1_9