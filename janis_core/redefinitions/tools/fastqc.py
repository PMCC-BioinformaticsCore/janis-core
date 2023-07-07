


from abc import ABC
from datetime import datetime
from typing import List, Dict, Any

from janis_core.operators.operator import IndexOperator
from janis_core import (
    ToolOutput,
    ToolInput,
    ToolMetadata,
    File,
    String,
    Boolean,
    Directory,
    Int,
    CaptureType,
    CpuSelector,
    InputSelector,
)
from janis_core import get_value_for_hints_and_ordered_resource_tuple
from janis_core.tool.test_classes import TTestCase
from ..types import ZipFile, TextFile, HtmlFile, FastqGz, FastqGzPair
from .bioinformaticstool import BioinformaticsTool


CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 2,
            CaptureType.EXOME: 2,
            CaptureType.CHROMOSOME: 4,
            CaptureType.THIRTYX: 8,
            CaptureType.NINETYX: 12,
            CaptureType.THREEHUNDREDX: 12,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 2,
            CaptureType.EXOME: 4,
            CaptureType.CHROMOSOME: 8,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 16,
            CaptureType.THREEHUNDREDX: 32,
        },
    )
]


class FastQCBase(BioinformaticsTool, ABC):
    def friendly_name(self) -> str:
        return "FastQC"

    def tool(self):
        return "fastqc"

    def base_command(self):
        return "fastqc"

    def tool_provider(self):
        return "FastQC"

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin", "Jiaan Yu"],
            dateCreated=datetime(2019, 3, 25),
            dateUpdated=datetime(2021, 11, 10),
            institution="Babraham Bioinformatics",
            doi=None,
            citation=None,
            keywords=["fastqc", "quality", "qa"],
            documentationUrl="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            documentation="FastQC is a program designed to spot potential problems in high througput sequencing datasets. "
            "It runs a set of analyses on one or more raw sequence files in fastq or bam format and produces a "
            "report which summarises the results.\n"
            "FastQC will highlight any areas where this library looks unusual and where you should take a closer look. "
            "The program is not tied to any specific type of sequencing technique and can be used to look at libraries "
            "coming from a large number of different experiment types "
            "(Genomic Sequencing, ChIP-Seq, RNA-Seq, BS-Seq etc etc).",
        )

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 1

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 8

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("reads", FastqGzPair),
            ToolInput(
                "read1",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 0),
                position=5,
            ),
            ToolInput(
                "read2",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 1),
                position=6,
            ),
            *self.additional_inputs,
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput(
                "out_R1",
                ZipFile(),
                selector=InputSelector("read1", remove_file_extension=True)
                + "_fastqc.zip",
            ),
            ToolOutput(
                "out_R1_datafile",
                TextFile,
                selector=InputSelector("read1", remove_file_extension=True)
                + "_fastqc/fastqc_data.txt",
            ),
            ToolOutput(
                "out_R1_html",
                HtmlFile,
                selector=InputSelector("read1", remove_file_extension=True)
                + "_fastqc.html",
            ),
            ToolOutput(
                "out_R1_directory",
                Directory,
                selector=InputSelector("read1", remove_file_extension=True) + "_fastqc",
            ),
            ToolOutput(
                "out_R2",
                ZipFile(),
                selector=InputSelector("read2", remove_file_extension=True)
                + "_fastqc.zip",
            ),
            ToolOutput(
                "out_R2_datafile",
                TextFile,
                selector=InputSelector("read2", remove_file_extension=True)
                + "_fastqc/fastqc_data.txt",
            ),
            ToolOutput(
                "out_R2_html",
                HtmlFile,
                selector=InputSelector("read2", remove_file_extension=True)
                + "_fastqc.html",
            ),
            ToolOutput(
                "out_R2_directory",
                Directory,
                selector=InputSelector("read2", remove_file_extension=True) + "_fastqc",
            ),
        ]

    additional_inputs = [
        ToolInput(
            "outdir",
            String(optional=True),
            default=".",
            prefix="--outdir",
            doc="(-o) Create all output files in the specified output directory. Please note that this "
            "directory must exist as the program will not create it.  If this option is not set then "
            "the output file for each sequence file is created in the same directory as the sequence "
            "file which was processed.",
        ),
        ToolInput(
            "casava",
            Boolean(optional=True),
            prefix="--casava",
            doc="Files come from raw casava output. Files in the same sample group "
            "(differing only by the group number) will be analysed as a set rather than individually. "
            "Sequences with the filter flag set in the header will be excluded from the analysis. "
            "Files must have the same names given to them by casava (including being gzipped and "
            "ending with .gz) otherwise they won't be grouped together correctly.",
        ),
        ToolInput(
            "nano",
            Boolean(optional=True),
            prefix="--nano",
            doc="Files come from naopore sequences and are in fast5 format. In this mode you can pass in "
            "directories to process and the program will take in all fast5 files within those "
            "directories and produce a single output file from the sequences found in all files.",
        ),
        ToolInput(
            "nofilter",
            Boolean(optional=True),
            prefix="--nofilter",
            doc="If running with --casava then don't remove read flagged by casava as poor quality when "
            "performing the QC analysis.",
        ),
        ToolInput(
            "extract",
            Boolean(optional=True),
            prefix="--extract",
            default=True,
            doc="If set then the zipped output file will be uncompressed in the same directory after it has "
            "been created.  By default this option will be set if fastqc is run in non-interactive mode.",
        ),
        ToolInput(
            "java",
            String(optional=True),
            prefix="--java",
            doc="(-j) Provides the full path to the java binary you want to use to launch fastqc. "
            "If not supplied then java is assumed to be in your path.",
        ),
        ToolInput(
            "noextract",
            Boolean(optional=True),
            prefix="--noextract",
            doc="Do not uncompress the output file after creating it.  You should set this option if you do"
            "not wish to uncompress the output when running in non-interactive mode. ",
        ),
        ToolInput(
            "nogroup",
            Boolean(optional=True),
            prefix="--nogroup",
            doc="Disable grouping of bases for reads >50bp. "
            "All reports will show data for every base in the read. "
            "WARNING: Using this option will cause fastqc to crash and burn if you use it on "
            "really long reads, and your plots may end up a ridiculous size. You have been warned! ",
        ),
        ToolInput(
            "format",
            String(optional=True),
            prefix="--format",
            doc="(-f) Bypasses the normal sequence file format detection and forces the program to use the "
            "specified format.  Valid formats are bam,sam,bam_mapped,sam_mapped and fastq ",
        ),
        ToolInput(
            "threads",
            Int(optional=True),
            prefix="--threads",
            default=CpuSelector(),
            doc="(-t) Specifies the number of files which can be processed simultaneously. "
            "Each thread will be allocated 250MB of memory so you shouldn't run more threads than your "
            "available memory will cope with, and not more than 6 threads on a 32 bit machine",
        ),
        ToolInput(
            "contaminants",
            File(optional=True),
            prefix="--contaminants",
            doc="(-c) Specifies a non-default file which contains the list of contaminants to screen "
            "overrepresented sequences against. The file must contain sets of named contaminants in "
            "the form name[tab]sequence.  Lines prefixed with a hash will be ignored.",
        ),
        ToolInput(
            "adapters",
            File(optional=True),
            prefix="--adapters",
            doc="(-a) Specifies a non-default file which contains the list of adapter sequences which will "
            "be explicity searched against the library. The file must contain sets of named adapters in "
            "the form name[tab]sequence. Lines prefixed with a hash will be ignored.",
        ),
        ToolInput(
            "limits",
            File(optional=True),
            prefix="--limits",
            doc="(-l) Specifies a non-default file which contains a set of criteria which will be used to "
            "determine the warn/error limits for the various modules.  This file can also be used to "
            "selectively  remove some modules from the output all together. "
            "The format needs to mirror the default limits.txt file found in the Configuration folder.",
        ),
        ToolInput(
            "kmers",
            Int(optional=True),
            prefix="--kmers",
            doc="(-k) Specifies the length of Kmer to look for in the Kmer content module. "
            "Specified Kmer length must be between 2 and 10. Default length is 7 if not specified. ",
        ),
        ToolInput(
            "quiet",
            Boolean(optional=True),
            prefix="--quiet",
            doc="(-q) Supress all progress messages on stdout and only report errors.",
        ),
        ToolInput(
            "dir",
            String(optional=True),
            prefix="--dir",
            doc="(-d) Selects a directory to be used for temporary files written when generating report images."
            "Defaults to system temp directory if not specified.",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "reads": [
                        f"{remote_dir}/NA12878-BRCA1_R1.fastq.gz",
                        f"{remote_dir}/NA12878-BRCA1_R2.fastq.gz",
                    ],
                    "threads": 1,
                },
                output=TextFile.basic_test(
                    "out_R1_datafile",
                    80000,
                    line_count=2739,
                    md5="8e23d29e0859ba547f0aa616ca395a8f",
                )
                + TextFile.basic_test(
                    "out_R2_datafile",
                    80000,
                    line_count=2751,
                    md5="a623b0a610ce51cef1137c7cb542f773",
                )
                + ZipFile.basic_test("out_R1", 400000)
                + ZipFile.basic_test("out_R2", 400000)
                + HtmlFile.basic_test("out_R1_html", 600000)
                + HtmlFile.basic_test("out_R2_html", 600000),
            ),
            TTestCase(
                name="minimal",
                input={
                    "reads": [
                        f"{remote_dir}/NA12878-BRCA1_R1.fastq.gz",
                        f"{remote_dir}/NA12878-BRCA1_R2.fastq.gz",
                    ],
                    "threads": 1,
                },
                output=self.minimal_test(),
            ),
        ]


class FastQCVersion_0_11_8:
    def version(self):
        return "v0.11.8"

    def container(self):
        return "quay.io/biocontainers/fastqc:0.11.8--2"



class FastQC_0_11_8(FastQCVersion_0_11_8, FastQCBase):
    pass