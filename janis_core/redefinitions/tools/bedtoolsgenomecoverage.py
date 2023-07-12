

from abc import ABC
from datetime import date

from janis_core import (
    ToolOutput,
    ToolInput,
    Boolean,
    Int,
    String,
    File,
    Float,
    InputSelector,
    Filename,
)
from janis_core.tool.test_classes import (
    TTestCase,
)

from .bioinformaticstool import BioinformaticsTool
from ..types import Bam, Bed, TextFile



class BedToolsToolBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "bedtools"

    def memory(self, hints):
        return 8


class BedToolsGenomeCoverageBedBase(BedToolsToolBase, ABC):
    def bind_metadata(self):

        self.metadata.contributors = ["Jiaan Yu"]
        self.metadata.dateUpdated = date(2020, 4, 1)
        self.metadata.dateCreated = date(2020, 4, 1)
        self.metadata.doi = None
        self.metadata.citation = None
        self.metadata.keywords = ["bedtools", "genomecov", "genomeCoverageBed"]
        self.metadata.documentationUrl = (
            "https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html"
        )
        self.metadata.documentation = """bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome. Note: 1. If using BED/GFF/VCF, the input (-i) file must be grouped by chromosome. A simple sort -k 1,1 in.bed > in.sorted.bed will suffice. Also, if using BED/GFF/VCF, one must provide a genome file via the -g argument. 2. If the input is in BAM (-ibam) format, the BAM file must be sorted by position. Using samtools sort aln.bam aln.sorted will suffice."""

    def tool(self):
        return "bedtoolsgenomeCoverageBed"

    def friendly_name(self):
        return "BEDTools: genomeCoverageBed"

    def base_command(self):
        return ["genomeCoverageBed"]

    def inputs(self):
        return [
            *self.additional_inputs,
            ToolInput(
                "inputBam",
                Bam(optional=True),
                prefix="-ibam",
                doc="Input bam file. Note: BAM _must_ be sorted by position. A 'samtools sort <BAM>' should suffice.",
            ),
            ToolInput(
                "inputBed",
                File(optional=True),
                prefix="-iBed",
                doc="Input bed file. Must be grouped by chromosome. A simple 'sort -k 1,1 <BED> > <BED>.sorted' will suffice.",
            ),
            ToolInput(
                "inputFile",
                File(optional=True),
                prefix="-i",
                doc="Input file, can be gff/vcf.",
            ),
            ToolInput(
                "genome",
                File(optional=True),
                prefix="-g",
                doc="Genome file. The genome file should tab delimited and structured as follows: <chromName><TAB><chromSize>.",
            ),
            ToolInput("outputFilename", Filename, prefix=">", position=10),
        ]

    def outputs(self):
        return [ToolOutput("out", TextFile, selector=InputSelector("outputFilename"))]

    additional_inputs = [
        ToolInput(
            "depth",
            Boolean(optional=True),
            prefix="-d",
            doc="Report the depth at each genome position (with one-based coordinates). Default behavior is to report a histogram.",
        ),
        ToolInput(
            "depthZero",
            Boolean(optional=True),
            prefix="-dz",
            doc="Report the depth at each genome position (with zero-based coordinates). Reports only non-zero positions. Default behavior is to report a histogram.",
        ),
        ToolInput(
            "BedGraphFormat",
            Boolean(optional=True),
            prefix="-bg",
            doc="Report depth in BedGraph format. For details, see: genome.ucsc.edu/goldenPath/help/bedgraph.html",
        ),
        ToolInput(
            "BedGraphFormata",
            Boolean(optional=True),
            prefix="-bga",
            doc="Report depth in BedGraph format, as above (-bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0  coverage by applying: 'grep -w 0$' to the output.",
        ),
        ToolInput(
            "split",
            Boolean(optional=True),
            prefix="-split",
            doc="Treat 'split' BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR 'N' and 'D' operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).",
        ),
        ToolInput(
            "strand",
            String(optional=True),
            prefix="-strand",
            doc="(STRING): can be + or -. Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).",
        ),
        ToolInput(
            "pairEnd",
            Boolean(optional=True),
            prefix="-pc",
            doc="Calculate coverage of pair-end fragments. Works for BAM files only",
        ),
        ToolInput(
            "fragmentSize",
            Boolean(optional=True),
            prefix="-fs",
            doc="Force to use provided fragment size instead of read length. Works for BAM files only",
        ),
        ToolInput(
            "du",
            Boolean(optional=True),
            prefix="-du",
            doc="Change strand af the mate read (so both reads from the same strand) useful for strand specific. Works for BAM files only",
        ),
        ToolInput(
            "fivePos",
            Boolean(optional=True),
            prefix="-5",
            doc="Calculate coverage of 5' positions (instead of entire interval).",
        ),
        ToolInput(
            "threePos",
            Boolean(optional=True),
            prefix="-3",
            doc="Calculate coverage of 3' positions (instead of entire interval).",
        ),
        ToolInput(
            "max",
            Int(optional=True),
            prefix="-max",
            doc="Combine all positions with a depth >= max into a single bin in the histogram. Irrelevant for -d and -bedGraph",
        ),
        ToolInput(
            "scale",
            Float(optional=True),
            prefix="-scale",
            doc="Scale the coverage by a constant factor. Each coverage value is multiplied by this factor before being reported. Useful for normalizing coverage by, e.g., reads per million (RPM). Default is 1.0; i.e., unscaled.",
        ),
        ToolInput(
            "trackline",
            Boolean(optional=True),
            prefix="-trackline",
            doc="Adds a UCSC/Genome-Browser track line definition in the first line of the output. - See here for more details about track line definition: http://genome.ucsc.edu/goldenPath/help/bedgraph.html - NOTE: When adding a trackline definition, the output BedGraph can be easily uploaded to the Genome Browser as a custom track, BUT CAN NOT be converted into a BigWig file (w/o removing the first line).",
        ),
        ToolInput(
            "trackopts",
            String(optional=True),
            prefix="-trackopts",
            doc="Writes additional track line definition parameters in the first line. - Example: -trackopts 'name=\"My Track\" visibility=2 color=255,30,30' Note the use of single-quotes if you have spaces in your parameters.",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "inputBam": f"{remote_dir}/NA12878-BRCA1.markduped.bam.bam",
                    "genome": f"{remote_dir}/NA12878-BRCA1.genome_file.txt",
                },
                output=TextFile.basic_test(
                    "out",
                    7432,
                    "chr17\t0\t83144233\t83257441\t0.99864",
                    220,
                    "f2007353bbd18f0a04eae9499d7c6a91",
                ),
            )
        ]



class BedTools_2_29_2:
    def container(self):
        return "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"

    def version(self):
        return "v2.29.2"


class BedToolsGenomeCoverageBed_2_29_2(BedTools_2_29_2, BedToolsGenomeCoverageBedBase):
    pass


BedToolsGenomeCoverageBedLatest = BedToolsGenomeCoverageBed_2_29_2