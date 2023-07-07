

from datetime import datetime
from typing import List

from janis_core import (
    ToolOutput,
    ToolInput,
    Boolean,
    String,
    Float,
    Int,
    Filename,
    ToolMetadata,
    get_value_for_hints_and_ordered_resource_tuple,
    InputSelector,
    CaptureType,
    Array,
)
from janis_core.tool.test_classes import TTestCase

from ..types import FastqGzPair
from ..tools import BioinformaticsTool


MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.EXOME: 4,
            CaptureType.CHROMOSOME: 8,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 16,
            CaptureType.THREEHUNDREDX: 24,
        },
    )
]

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 2,
            CaptureType.EXOME: 4,
            CaptureType.CHROMOSOME: 5,
            CaptureType.THIRTYX: 8,
            CaptureType.NINETYX: 12,
            CaptureType.THREEHUNDREDX: 16,
        },
    )
]


class CutAdaptBase_2(BioinformaticsTool):
    def friendly_name(self) -> str:
        return "Cutadapt"

    def tool_provider(self):
        return "cutadapt"

    def tool(self) -> str:
        return "cutadapt"

    def base_command(self):
        return "cutadapt"

    def inputs(self) -> List[ToolInput]:
        import uuid

        return [
            ToolInput("outputPrefix", String(), doc="Used for naming purposes"),
            ToolInput("fastq", FastqGzPair, position=5),
            ToolInput(
                "adapter",
                input_type=Array(String(), optional=True),
                prefix="-a",
                prefix_applies_to_all_elements=True,
                doc="Sequence of an adapter ligated to the 3' end (paired data: of the first read). "
                "The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), "
                "the adapter is only found if it is a suffix of the read.",
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("outputPrefix"),
                    suffix="R1",
                    extension=".fastq.gz",
                ),
                prefix="-o",
                doc="Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. "
                "The summary report is sent to standard output. Use '{name}' in FILE to demultiplex "
                "reads into multiple files. Default: write to standard output",
            ),
            ToolInput(
                "secondReadFile",
                Filename(
                    prefix=InputSelector("outputPrefix"),
                    suffix="R2",
                    extension=".fastq.gz",
                ),
                prefix="-p",
                doc="Write second read in a pair to FILE.",
            ),
            *self.additional_args,
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput(
                "out",
                FastqGzPair,
                selector=[
                    InputSelector("outputPrefix") + "-R1.fastq.gz",
                    InputSelector("outputPrefix") + "-R2.fastq.gz",
                ],
            )
        ]

    def memory(self, hints):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 4

    def cpus(self, hints):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 5

    additional_args = [
        ToolInput(
            tag="cores",
            input_type=Int(optional=True),
            prefix="--cores",
            separate_value_from_prefix=True,
            doc="(-j)  Number of CPU cores to use. Use 0 to auto-detect. Default: 1",
        ),
        ToolInput(
            tag="front",
            input_type=Array(String(), optional=True),
            prefix="--front",
            separate_value_from_prefix=True,
            doc="(-g)  Sequence of an adapter ligated to the 5' end (paired data: of the first read). The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only found if it is a prefix of the read.",
        ),
        ToolInput(
            tag="anywhere",
            input_type=Array(String(), optional=True),
            prefix="--anywhere",
            separate_value_from_prefix=True,
            doc="(-b)  Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of the first read). Both types of matches as described under -a und -g are allowed. If the first base of the read is part of the match, the behavior is as with -g, otherwise as with -a. This option is mostly for rescuing failed library preparations - do not use if you know which end your adapter was ligated to!",
        ),
        ToolInput(
            tag="errorRate",
            input_type=Float(optional=True),
            prefix="--error-rate",
            separate_value_from_prefix=True,
            doc="(-e)  Maximum allowed error rate as value between 0 and 1 (no. of errors divided by length of matching region). Default: 0.1 (=10%)",
        ),
        ToolInput(
            tag="noIndels",
            input_type=Boolean(optional=True),
            prefix="--no-indels",
            separate_value_from_prefix=True,
            doc="Allow only mismatches in alignments. Default: allow both mismatches and indels",
        ),
        ToolInput(
            tag="times",
            input_type=Int(optional=True),
            prefix="--times",
            separate_value_from_prefix=True,
            doc="(-n)  Remove up to COUNT adapters from each read. Default: 1",
        ),
        ToolInput(
            tag="overlap",
            input_type=Int(optional=True),
            prefix="--overlap",
            separate_value_from_prefix=True,
            doc="(-O)  Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3",
        ),
        ToolInput(
            tag="matchReadWildcards",
            input_type=Boolean(optional=True),
            prefix="--match-read-wildcards",
            separate_value_from_prefix=True,
            doc=" Interpret IUPAC wildcards in reads. Default: False",
        ),
        ToolInput(
            tag="noMatchAdapterWildcards",
            input_type=Boolean(optional=True),
            prefix="--no-match-adapter-wildcards",
            separate_value_from_prefix=True,
            doc="(-N)  Do not interpret IUPAC wildcards in adapters.",
        ),
        ToolInput(
            tag="action",
            input_type=String(optional=True),
            prefix="--action",
            separate_value_from_prefix=True,
            doc="(trim,mask,lowercase,none}  What to do with found adapters. mask: replace with 'N' characters; lowercase: convert to lowercase; none: leave unchanged (useful with --discard-untrimmed). Default: trim",
        ),
        ToolInput(
            tag="cut",
            input_type=Int(optional=True),
            prefix="--cut",
            separate_value_from_prefix=True,
            doc="(-u)  Remove bases from each read (first read only if paired). If LENGTH is positive, remove bases from the beginning. If LENGTH is negative, remove bases from the end. Can be used twice if LENGTHs have different signs. This is applied *before* adapter trimming.",
        ),
        ToolInput(
            tag="nextseqTrim",
            input_type=String(optional=True),
            prefix="--nextseq-trim",
            separate_value_from_prefix=True,
            doc=" NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.",
        ),
        ToolInput(
            tag="qualityCutoff",
            input_type=Int(optional=True),
            prefix="--quality-cutoff",
            separate_value_from_prefix=True,
            doc="(]3'CUTOFF, ]3'CUTOFF, -q)  Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.",
        ),
        ToolInput(
            tag="qualityBase",
            input_type=Boolean(optional=True),
            prefix="--quality-base",
            separate_value_from_prefix=True,
            doc="Assume that quality values in FASTQ are encoded as ascii(quality + N). This needs to be set to 64 for some old Illumina FASTQ files. Default: 33",
        ),
        ToolInput(
            tag="length",
            input_type=Int(optional=True),
            prefix="--length",
            separate_value_from_prefix=True,
            doc="(-l)  Shorten reads to LENGTH. Positive values remove bases at the end while negative ones remove bases at the beginning. This and the following modifications are applied after adapter trimming.",
        ),
        ToolInput(
            tag="trimN",
            input_type=Int(optional=True),
            prefix="--trim-n",
            separate_value_from_prefix=True,
            doc="Trim N's on ends of reads.",
        ),
        ToolInput(
            tag="lengthTag",
            input_type=Int(optional=True),
            prefix="--length-tag",
            separate_value_from_prefix=True,
            doc="Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.",
        ),
        ToolInput(
            tag="stripSuffix",
            input_type=String(optional=True),
            prefix="--strip-suffix",
            separate_value_from_prefix=True,
            doc=" Remove this suffix from read names if present. Can be given multiple times.",
        ),
        ToolInput(
            tag="prefix",
            input_type=String(optional=True),
            prefix="--prefix",
            separate_value_from_prefix=True,
            doc="(-x)  Add this prefix to read names. Use {name} to insert the name of the matching adapter.",
        ),
        ToolInput(
            tag="suffix",
            input_type=String(optional=True),
            prefix="--suffix",
            separate_value_from_prefix=True,
            doc="(-y)  Add this suffix to read names; can also include {name}",
        ),
        ToolInput(
            tag="zeroCap",
            input_type=Boolean(optional=True),
            prefix="--zero-cap",
            separate_value_from_prefix=True,
            doc="(-z) Change negative quality values to zero.",
        ),
        ToolInput(
            tag="minimumLength",
            input_type=Int(optional=True),
            prefix="--minimum-length",
            separate_value_from_prefix=True,
            doc="(-m)  Discard reads shorter than LEN. Default: 0",
        ),
        ToolInput(
            tag="maximumLength",
            input_type=Int(optional=True),
            prefix="--maximum-length",
            separate_value_from_prefix=True,
            doc="(-M)  Discard reads longer than LEN. Default: no limit",
        ),
        ToolInput(
            tag="maxN",
            input_type=Float(optional=True),
            prefix="--max-n",
            separate_value_from_prefix=True,
            doc="Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.",
        ),
        ToolInput(
            tag="discardTrimmed",
            input_type=Boolean(optional=True),
            prefix="--discard-trimmed",
            separate_value_from_prefix=True,
            doc="(--discard)  Discard reads that contain an adapter. Use also -O to avoid discarding too many randomly matching reads.",
        ),
        ToolInput(
            tag="discardUntrimmed",
            input_type=Boolean(optional=True),
            prefix="--discard-untrimmed",
            separate_value_from_prefix=True,
            doc="(--trimmed-only)  Discard reads that do not contain an adapter.",
        ),
        ToolInput(
            tag="discardCasava",
            input_type=Boolean(optional=True),
            prefix="--discard-casava",
            separate_value_from_prefix=True,
            doc="Discard reads that did not pass CASAVA filtering (header has :Y:).",
        ),
        ToolInput(
            tag="quiet",
            input_type=Boolean(optional=True),
            prefix="--quiet",
            separate_value_from_prefix=True,
            doc="Print only error messages. Which type of report to print: 'full' or 'minimal'. Default: full",
        ),
        ToolInput(
            tag="compressionLevel",
            input_type=String(optional=True),
            prefix="-Z",
            separate_value_from_prefix=True,
            doc="Use compression level 1 for gzipped output files (faster, but uses more space)",
        ),
        ToolInput(
            tag="infoFile",
            input_type=String(optional=True),
            prefix="--info-file",
            separate_value_from_prefix=True,
            doc="Write information about each read and its adapter matches into FILE. See the documentation for the file format.",
        ),
        ToolInput(
            tag="restFile",
            input_type=String(optional=True),
            prefix="--rest-file",
            separate_value_from_prefix=True,
            doc="(-r)  When the adapter matches in the middle of a read, write the rest (after the adapter) to FILE.",
        ),
        ToolInput(
            tag="wildcardFile",
            input_type=String(optional=True),
            prefix="--wildcard-file",
            separate_value_from_prefix=True,
            doc="When the adapter has N wildcard bases, write adapter bases matching wildcard positions to FILE. (Inaccurate with indels.)",
        ),
        ToolInput(
            tag="tooShortOutput",
            input_type=String(optional=True),
            prefix="--too-short-output",
            separate_value_from_prefix=True,
            doc=" Write reads that are too short (according to length specified by -m) to FILE. Default: discard reads",
        ),
        ToolInput(
            tag="tooLongOutput",
            input_type=String(optional=True),
            prefix="--too-long-output",
            separate_value_from_prefix=True,
            doc=" Write reads that are too long (according to length specified by -M) to FILE. Default: discard reads",
        ),
        ToolInput(
            tag="untrimmedOutput",
            input_type=String(optional=True),
            prefix="--untrimmed-output",
            separate_value_from_prefix=True,
            doc=" Write reads that do not contain any adapter to FILE. Default: output to same file as trimmed reads",
        ),
        ToolInput(
            tag="adapterSecondRead",
            input_type=Array(String, optional=True),
            prefix="-A",
            separate_value_from_prefix=True,
            prefix_applies_to_all_elements=True,
            doc="3' adapter to be removed from second read in a pair.",
        ),
        ToolInput(
            tag="frontAdapterSecondRead",
            input_type=Array(String, optional=True),
            prefix="-G",
            separate_value_from_prefix=True,
            doc="5' adapter to be removed from second read in a pair.",
        ),
        ToolInput(
            tag="anywhereAdapterSecondRead",
            input_type=Array(String, optional=True),
            prefix="-B",
            separate_value_from_prefix=True,
            doc="5'/3 adapter to be removed from second read in a pair.",
        ),
        ToolInput(
            tag="removeNBasesFromSecondRead",
            input_type=String(optional=True),
            prefix="-U",
            separate_value_from_prefix=True,
            doc="Remove LENGTH bases from second read in a pair.",
        ),
        ToolInput(
            tag="pairAdapters",
            input_type=Boolean(optional=True),
            prefix="--pair-adapters",
            separate_value_from_prefix=True,
            doc="Treat adapters given with -a/-A etc. as pairs. Either both or none are removed from each read pair.",
        ),
        ToolInput(
            tag="pairFilter",
            input_type=String(optional=True),
            prefix="--pair-filter",
            separate_value_from_prefix=True,
            doc="{any,both,first} Which of the reads in a paired-end read have to match the filtering criterion in order for the pair to be filtered. Default: any",
        ),
        ToolInput(
            tag="interleaved",
            input_type=Boolean(optional=True),
            prefix="--interleaved",
            separate_value_from_prefix=True,
            doc="Read and write interleaved paired-end reads.",
        ),
        ToolInput(
            tag="untrimmedPairedOutput",
            input_type=String(optional=True),
            prefix="--untrimmed-paired-output",
            separate_value_from_prefix=True,
            doc=" Write second read in a pair to this FILE when no adapter was found. Use with --untrimmed-output. Default: output to same file as trimmed reads",
        ),
        ToolInput(
            tag="tooShortPairedOutput",
            input_type=String(optional=True),
            prefix="--too-short-paired-output",
            separate_value_from_prefix=True,
            doc=" Write second read in a pair to this file if pair is too short. Use also --too-short-output.",
        ),
        ToolInput(
            tag="tooLongPairedOutput",
            input_type=String(optional=True),
            prefix="--too-long-paired-output",
            separate_value_from_prefix=True,
            doc=" Write second read in a pair to this file if pair is too long. Use also --too-long-output.",
        ),
    ]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin", "Jiaan Yu"],
            dateCreated=datetime(2019, 3, 21),
            dateUpdated=datetime(2021, 11, 3),
            institution=None,
            doi="DOI:10.14806/ej.17.1.200",
            citation="Martin, Marcel. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” "
            "EMBnet.journal, vol. 17, no. 1, EMBnet Stichting, May 2011, p. 10. "
            "Crossref, doi:10.14806/ej.17.1.200.",
            keywords=["cutadapt", "trim"],
            documentationUrl="https://cutadapt.readthedocs.io/en/stable/",
            documentation='cutadapt version 2.4\nCopyright (C) 2010-2019 Marcel Martin <marcel.martin@scilifelab.se>\ncutadapt removes adapter sequences from high-throughput sequencing reads.\nUsage:\n    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq\nFor paired-end reads:\n    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq\nReplace "ADAPTER" with the actual sequence of your 3\' adapter. IUPAC wildcard\ncharacters are supported. The reverse complement is *not* automatically\nsearched. All reads from input.fastq will be written to output.fastq with the\nadapter sequence removed. Adapter matching is error-tolerant. Multiple adapter\nsequences can be given (use further -a options), but only the best-matching\nadapter will be removed.\nInput may also be in FASTA format. Compressed input and output is supported and\nauto-detected from the file name (.gz, .xz, .bz2). Use the file name \'-\' for\nstandard input/output. Without the -o option, output is sent to standard output.\nCitation:\nMarcel Martin. Cutadapt removes adapter sequences from high-throughput\nsequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.\nhttp://dx.doi.org/10.14806/ej.17.1.200\nRun "cutadapt - -help" to see all command-line options.\nSee https://cutadapt.readthedocs.io/ for full documentation.\n',
        )

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "fastq": [
                        f"{remote_dir}/NA12878-BRCA1_R1.fastq.gz",
                        f"{remote_dir}/NA12878-BRCA1_R2.fastq.gz",
                    ],
                    "qualityCutoff": 15,
                    "minimumLength": 50,
                    "outputPrefix": "output",
                },
                output=FastqGzPair.basic_test(
                    "out",
                    1090240,
                    1163374,
                ),
            ),
            TTestCase(
                name="minimal",
                input={
                    "fastq": [
                        f"{remote_dir}/NA12878-BRCA1_R1.fastq.gz",
                        f"{remote_dir}/NA12878-BRCA1_R2.fastq.gz",
                    ],
                    "qualityCutoff": 15,
                    "minimumLength": 50,
                    "outputPrefix": "output",
                },
                output=self.minimal_test(),
            ),
        ]


class CutAdapt_2_1(CutAdaptBase_2):
    def container(self):
        return "quay.io/biocontainers/cutadapt:2.1--py37h14c3975_0"

    def version(self):
        return "2.1"


class CutAdapt_2_4(CutAdaptBase_2):
    def container(self):
        return "quay.io/biocontainers/cutadapt:2.4--py37h14c3975_0"

    def version(self):
        return "2.4"


class CutAdapt_2_5(CutAdaptBase_2):
    def container(self):
        return "quay.io/biocontainers/cutadapt:2.5--py37h516909a_0"

    def version(self):
        return "2.5"


class CutAdapt_2_6(CutAdaptBase_2):
    def container(self):
        return "quay.io/biocontainers/cutadapt:2.6--py36h516909a_0"

    def version(self):
        return "2.6"


CutAdaptLatest = CutAdapt_2_6