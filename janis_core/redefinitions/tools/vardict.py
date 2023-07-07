
from datetime import datetime
from abc import ABC
from typing import List, Dict, Any

from janis_core import (
    CpuSelector,
    TOutput,
    ToolOutput,
    ToolInput,
    Filename,
    File,
    OutputDocumentation,
    ToolArgument,
    Boolean,
    Float,
    Int,
    String,
    InputSelector,
    CaptureType,
    CommandTool,
    get_value_for_hints_and_ordered_resource_tuple,
    ToolMetadata,
)

from ..types import FastaDict, BamBai, Bed, FastaFai, Vcf
from .bioinformaticstool import BioinformaticsPythonTool, BioinformaticsTool


### VERSIONS ###



class VarDict_1_5_6(CommandTool, ABC):
    def container(self):
        return "michaelfranklin/vardict:1.5.6"

    def version(self):
        return "1.5.6"


class VarDict_1_5_7(CommandTool):
    def container(self):
        return "michaelfranklin/vardict:1.5.7"

    def version(self):
        return "1.5.7"


class VarDict_1_5_8(CommandTool):
    def container(self):
        return "michaelfranklin/vardict:1.5.8"

    def version(self):
        return "1.5.8"


class VarDict_1_6_0(CommandTool):
    def container(self):
        return "michaelfranklin/vardict:1.6.0"

    def version(self):
        return "1.6.0"


class VarDict_1_7_0(CommandTool):
    def container(self):
        return "michaelfranklin/vardict:1.7.0"

    def version(self):
        return "1.7.0"



### GENERATE HEADER LINES ###


class GenerateVardictHeaderLines(BioinformaticsPythonTool):
    @staticmethod
    def code_block(
        reference: FastaDict, output_filename: str = "output.txt"
    ) -> Dict[str, Any]:
        """
        :param reference: Reference file to generate vardict header lines for (must have ^.dict) pattern
        :param output_filename: Filename to output to
        """
        from re import sub

        ref_dict = sub("\.fa(sta)?$", ".dict", reference)

        with open(output_filename, "w+") as out, open(ref_dict) as inp:
            out.write("##source=vardict\n")
            for line in inp:
                if not line.startswith("@SQ"):
                    continue
                pieces = line.split("\t")
                chrom = pieces[1].replace("SN:", "")
                length = pieces[2].replace("LN:", "")

                out.write(f"##contig=<ID={chrom},length={length}>\n")

            return {"out": output_filename}

    def outputs(self) -> List[TOutput]:
        return [
            TOutput(
                "out",
                File,
                doc=OutputDocumentation(
                    doc="Header file for VarDict, generated based on the reference index"
                ),
            )
        ]

    def id(self) -> str:
        return "GenerateVardictHeaderLines"

    def friendly_name(self) -> str:
        return "GenerateVardictHeaderLines"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def version(self):
        return "v0.1.0"

    def bind_metadata(self):
        self.metadata.dateCreated = datetime(2020, 6, 2)
        self.metadata.dateUpdated = datetime(2020, 6, 2)
        self.metadata.contributors = ["Michael Franklin", "Jiaan Yu"]
        self.metadata.documentation = """\
Generate VarDict Headerlines.       
        """




### VARDICT GERMLINE ###


CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.CHROMOSOME: 8,
            CaptureType.EXOME: 8,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 16,
            CaptureType.THREEHUNDREDX: 16,
        },
    )
]


MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 32,
            CaptureType.THIRTYX: 64,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class VarDictGermlineBase(BioinformaticsTool, ABC):
    def friendly_name(self) -> str:
        return "VarDict (Germline)"

    def tool(self):
        return "vardict_germline"

    def tool_provider(self):
        return "VarDict"

    def base_command(self):
        return "VarDict"

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 4

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 8

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("intervals", Bed(), position=2, shell_quote=False),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf", suffix=".vardict"),
                prefix=">",
                position=6,
                shell_quote=False,
            ),
            ToolInput(
                "bam",
                BamBai(),
                prefix="-b",
                position=1,
                shell_quote=False,
                doc="The indexed BAM file",
            ),
            ToolInput(
                "reference",
                FastaFai(),
                prefix="-G",
                position=1,
                shell_quote=False,
                doc="The reference fasta. Should be indexed (.fai). "
                "Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa",
            ),
            *VarDictGermlineBase.vardict_inputs,
            *VarDictGermlineBase.var2vcf_inputs,
        ]

    def outputs(self):
        return [ToolOutput("out", Vcf(), glob=InputSelector("outputFilename"))]

    def arguments(self):
        return [
            # ToolArgument("export VarDict=\"/config/binaries/vardict/1.5.1/bin/VarDict\";", position=0, shell_quote=False),
            # ToolArgument("", position=0, shell_quote=False),
            ToolArgument("| teststrandbias.R |", position=3, shell_quote=False),
            ToolArgument("var2vcf_valid.pl", position=4, shell_quote=False),
        ]

    vardict_inputs = [
        ToolInput(
            "indels3prime",
            Boolean(optional=True),
            prefix="-3",
            position=1,
            shell_quote=False,
            doc="Indicate to move indels to 3-prime if alternative alignment can be achieved.",
        ),
        ToolInput(
            "amplicon",
            Float(optional=True),
            prefix="-a",
            position=1,
            shell_quote=False,
            doc="Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  "
            "A read pair is considered belonging  to the amplicon if the edges are less than int bp to "
            "the amplicon, and overlap fraction is at least float.  Default: 10:0.95",
        ),
        ToolInput(
            "minReads",
            Int(optional=True),
            prefix="-B",
            position=1,
            shell_quote=False,
            doc="The minimum # of reads to determine strand bias, default 2",
        ),
        ToolInput(
            "chromNamesAreNumbers",
            Boolean(optional=True),
            prefix="-C",
            position=1,
            shell_quote=False,
            doc="Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2",
        ),
        ToolInput(
            "chromColumn",
            Int(optional=True),
            prefix="-c",
            position=1,
            shell_quote=False,
            doc="The column for chromosome",
        ),
        ToolInput(
            "debug",
            Boolean(optional=True),
            prefix="-D",
            position=1,
            shell_quote=False,
            doc="Debug mode.  Will print some error messages and append full genotype at the end.",
        ),
        ToolInput(
            "splitDelimeter",
            String(optional=True),
            prefix="-d",
            position=1,
            shell_quote=False,
            doc='The delimiter for split region_info, default to tab "\t"',
        ),
        ToolInput(
            "geneEndCol",
            Int(optional=True),
            prefix="-E",
            position=1,
            shell_quote=False,
            doc="The column for region end, e.g. gene end",
        ),
        ToolInput(
            "segEndCol",
            Int(optional=True),
            prefix="-e",
            position=1,
            shell_quote=False,
            doc="The column for segment ends in the region, e.g. exon ends",
        ),
        ToolInput(
            "filter",
            String(optional=True),
            prefix="-F",
            position=1,
            shell_quote=False,
            doc="The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and "
            "duplicates). Use -F 0 to turn it off.",
        ),
        ToolInput(
            "alleleFreqThreshold",
            Float(optional=True),
            prefix="-f",
            position=1,
            shell_quote=False,
            doc="The threshold for allele frequency, default: 0.05 or 5%",
        ),
        ToolInput(
            "geneNameCol",
            Int(optional=True),
            prefix="-g",
            position=1,
            shell_quote=False,
            doc="The column for gene name, or segment annotation",
        ),
        # ToolInput("help", Boolean(optional=True), prefix="-H", position=1, shell_quote=False,
        #           doc="Print this help page"),
        ToolInput(
            "printHeaderRow",
            Boolean(optional=True),
            prefix="-h",
            position=1,
            shell_quote=False,
            doc="Print a header row describing columns",
        ),
        ToolInput(
            "indelSize",
            Int(optional=True),
            prefix="-I",
            position=1,
            shell_quote=False,
            doc="The indel size.  Default: 120bp",
        ),
        ToolInput(
            "outputSplice",
            Boolean(optional=True),
            prefix="-i",
            position=1,
            shell_quote=False,
            doc="Output splicing read counts",
        ),
        ToolInput(
            "performLocalRealignment",
            Int(optional=True),
            prefix="-k",
            position=1,
            shell_quote=False,
            doc="Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. "
            "For Ion or PacBio, 0 is recommended.",
        ),
        ToolInput(
            "minMatches",
            Int(optional=True),
            prefix="-M",
            position=1,
            shell_quote=False,
            doc="The minimum matches for a read to be considered. If, after soft-clipping, the matched "
            "bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing "
            "where there's no insert and the matching is only the primers. Default: 0, or no filtering",
        ),
        ToolInput(
            "maxMismatches",
            Int(optional=True),
            prefix="-m",
            position=1,
            shell_quote=False,
            doc="If set, reads with mismatches more than INT will be filtered and ignored. "
            "Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln "
            "followed by sampe. BWA mem is calculated as NM - Indels. "
            "Default: 8, or reads with more than 8 mismatches will not be used.",
        ),
        ToolInput(
            "sampleName",
            String(),
            prefix="-N",
            position=1,
            shell_quote=False,
            doc="The sample name to be used directly.  Will overwrite -n option",
        ),
        ToolInput(
            "regexSampleName",
            String(optional=True),
            prefix="-n",
            position=1,
            shell_quote=False,
            doc="The regular expression to extract sample name from BAM filenames. "
            "Default to: /([^\/\._]+?)_[^\/]*.bam/",
        ),
        ToolInput(
            "mapq",
            String(optional=True),
            prefix="-O",
            position=1,
            shell_quote=False,
            doc="The reads should have at least mean MapQ to be considered a valid variant. "
            "Default: no filtering",
        ),
        ToolInput(
            "qratio",
            Float(optional=True),
            prefix="-o",
            position=1,
            shell_quote=False,
            doc="The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). "
            "The quality is defined by -q option.  Default: 1.5",
        ),
        ToolInput(
            "readPosition",
            Float(optional=True),
            prefix="-P",
            position=1,
            shell_quote=False,
            doc="The read position filter. If the mean variants position is less that specified, "
            "it's considered false positive.  Default: 5",
        ),
        ToolInput(
            "pileup",
            Boolean(optional=True),
            prefix="-p",
            position=1,
            shell_quote=False,
            doc="Do pileup regardless of the frequency",
        ),
        ToolInput(
            "minMappingQual",
            Int(optional=True),
            prefix="-Q",
            position=1,
            shell_quote=False,
            doc="If set, reads with mapping quality less than INT will be filtered and ignored",
        ),
        ToolInput(
            "phredScore",
            Int(optional=True),
            prefix="-q",
            position=1,
            shell_quote=False,
            doc="The phred score for a base to be considered a good call.  "
            "Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.",
        ),
        ToolInput(
            "region",
            String(optional=True),
            prefix="-R",
            position=1,
            shell_quote=False,
            doc="The region of interest.  In the format of chr:start-end.  If end is omitted, "
            "then a single position.  No BED is needed.",
        ),
        ToolInput(
            "minVariantReads",
            Int(optional=True),
            prefix="-r",
            position=1,
            shell_quote=False,
            doc="The minimum # of variant reads, default 2",
        ),
        ToolInput(
            "regStartCol",
            Int(optional=True),
            prefix="-S",
            position=1,
            shell_quote=False,
            doc="The column for region start, e.g. gene start",
        ),
        ToolInput(
            "segStartCol",
            Int(optional=True),
            prefix="-s",
            position=1,
            shell_quote=False,
            doc="The column for segment starts in the region, e.g. exon starts",
        ),
        ToolInput(
            "minReadsBeforeTrim",
            Int(optional=True),
            prefix="-T",
            position=1,
            shell_quote=False,
            doc="Trim bases after [INT] bases in the reads",
        ),
        ToolInput(
            "removeDuplicateReads",
            Boolean(optional=True),
            prefix="-t",
            position=1,
            shell_quote=False,
            doc="Indicate to remove duplicated reads.  Only one pair with same start positions will be kept",
        ),
        ToolInput(
            "threads",
            Int(optional=True),
            default=CpuSelector(),
            prefix="-th",
            position=1,
            shell_quote=False,
            doc="Threads count.",
        ),
        ToolInput(
            "freq",
            Int(optional=True),
            prefix="-V",
            position=1,
            shell_quote=False,
            doc="The lowest frequency in the normal sample allowed for a putative somatic mutation. "
            "Defaults to 0.05",
        ),
        ToolInput(
            "vcfFormat",
            Boolean(optional=True),
            prefix="-v",
            position=1,
            shell_quote=False,
            doc="VCF format output",
        ),
        ToolInput(
            "vs",
            String(optional=True),
            prefix="-VS",
            position=1,
            shell_quote=False,
            doc="[STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: "
            "STRICT   - throw an exception if something looks wrong. "
            "LENIENT	- Emit warnings but keep going if possible. "
            "SILENT	- Like LENIENT, only don't emit warning messages. "
            "Default: LENIENT",
        ),
        ToolInput(
            "bp",
            Int(optional=True),
            prefix="-X",
            position=1,
            shell_quote=False,
            doc="Extension of bp to look for mismatches after insersion or deletion.  "
            "Default to 3 bp, or only calls when they're within 3 bp.",
        ),
        ToolInput(
            "extensionNucleotide",
            Int(optional=True),
            prefix="-x",
            position=1,
            shell_quote=False,
            doc="The number of nucleotide to extend for each segment, default: 0",
        ),
        ToolInput(
            "yy",
            Boolean(optional=True),
            prefix="-y",
            position=1,
            shell_quote=False,
            doc="<No content>",
        ),
        ToolInput(
            "downsamplingFraction",
            Int(optional=True),
            prefix="-Z",
            position=1,
            shell_quote=False,
            doc="For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  "
            "Default: No downsampling.  Use with caution.  "
            "The downsampling will be random and non-reproducible.",
        ),
        ToolInput(
            "zeroBasedCoords",
            Int(optional=True),
            prefix="-z",
            position=1,
            shell_quote=False,
            doc="0/1  Indicate whether coordinates are zero-based, as IGV uses.  "
            "Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. "
            "When using the -R option, it's set to 0",
        ),
    ]

    var2vcf_inputs = [
        ToolInput(
            "var2vcfSampleName", String(), prefix="-N", position=5, shell_quote=False
        ),
        ToolInput(
            "var2vcfAlleleFreqThreshold",
            Float(),
            prefix="-f",
            position=5,
            shell_quote=False,
        ),
    ]

    def docurl(self):
        return "https://github.com/AstraZeneca-NGS/VarDict"

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 1, 21),
            dateUpdated=datetime(2020, 6, 4),
            documentation="",
        )

    def doc(self):
        return """
    VarDict
    
    VarDict is an ultra sensitive variant caller for both single and paired sample variant 
    calling from BAM files. VarDict implements several novel features such as amplicon bias 
    aware variant calling from targeted sequencing experiments, rescue of long indels by 
    realigning bwa soft clipped reads and better scalability than many Java based variant callers.

    Due to the philosophy of VarDict in calling "everything", several downstream strategies have 
    been developed to filter variants to for example the most likely cancer driving events. 
    These strategies are based on evidence in different databases and/or quality metrics. 
    http://bcb.io/2016/04/04/vardict-filtering/ provides an overview of how to develop further 
    filters for VarDict. The script at https://github.com/AstraZeneca-NGS/VarDict/blob/master/vcf2txt.pl 
    can be used to put the variants into a context by including information from dbSNP, Cosmic and ClinVar. 
    We are open to suggestions from the community on how to best narrow down to the variants of most interest.
    
    A Java based drop-in replacement for vardict.pl is being developed at 
    https://github.com/AstraZeneca-NGS/VarDictJava. The Java implementation is approximately 
    10 times faster than the original Perl implementation and does not depend on samtools
    
    To enable amplicon aware variant calling (single sample mode only; not supported in paired 
    variant calling), please make sure the bed file has 8 columns with the 7th and 8th columns 
    containing the insert interval (therefore subset of the 2nd and 3rd column interval). 
    
    Requirements

        - Perl (uses /usr/bin/env perl)
        - R (uses /usr/bin/env R)
        - samtools (must be in path, not required if using the Java implementation in place of vardict.pl)
    """


class VarDictGermline_1_5_6(VarDictGermlineBase, VarDict_1_5_6):
    pass


class VarDictGermline_1_5_7(VarDictGermlineBase, VarDict_1_5_7):
    pass


class VarDictGermline_1_5_8(VarDictGermlineBase, VarDict_1_5_8):
    pass


class VarDictGermline_1_6_0(VarDictGermlineBase, VarDict_1_6_0):
    pass


class VarDictGermline_1_7_0(VarDictGermlineBase, VarDict_1_7_0):
    pass
