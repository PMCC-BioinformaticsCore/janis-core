from datetime import date
from typing import Dict, Any

from janis_core import (
    ToolOutput,
    ToolInput,
    Filename,
    File,
    String,
    Float,
    Int,
    Boolean,
    Array,
    InputSelector,
    CpuSelector,
    CommandTool
)

from janis_core.redefinitions.types import Bam, BamBai, FastaWithDict, Bed, Vcf


class GridssTestTool(CommandTool):
    def tool(self) -> str:
        return "gridss"
    
    def container(self):
        return "gridss/gridss:2.10.2"

    def version(self) -> str:
        return "v2.10.2"

    def tool_provider(self):
        return "Papenfuss Labs"

    def friendly_name(self) -> str:
        return "Gridss"

    def base_command(self):
        return "/opt/gridss/gridss.sh"

    def inputs(self):
        return [
            ToolInput("bams", Array(BamBai()), position=10),
            ToolInput(
                "reference",
                FastaWithDict(),
                prefix="--reference",
                doc="reference genome to use.",
            ),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf"),
                prefix="--output",
                doc="output VCF.",
            ),
            ToolInput(
                "assemblyFilename",
                Filename(suffix=".assembly", extension=".bam"),
                prefix="--assembly",
                doc="location of the GRIDSS assembly BAM. This file will be created by GRIDSS.",
            ),
            ToolInput(
                "threads",
                Int(optional=True),
                default=CpuSelector(),
                prefix="--threads",
                doc="number of threads to use. (Default: 8)",
            ),
            ToolInput(
                "jarPath",
                String(optional=True),
                prefix="--jar",
                doc="location of GRIDSS jar",
            ),
            ToolInput(
                "workingDir",
                String(optional=True),
                default="./TMP",
                prefix="--workingdir",
                doc="directory to place GRIDSS intermediate and temporary files. .gridss.working subdirectories will be created. (Default: .)",
            ),
            ToolInput(
                "blacklist",
                Bed(optional=True),
                prefix="--blacklist",
                doc="BED file containing regions to ignore",
            ),
            ToolInput(
                "steps",
                Array(String, optional=True),
                prefix="--steps",
                separator=",",
                prefix_applies_to_all_elements=False,
                doc="processing steps to run. Defaults to all steps. Multiple steps are specified using comma separators. Possible steps are: setupreference, preprocess, assemble, call, all. WARNING: multiple instances of GRIDSS generating reference files at the same time will result in file corruption. Make sure these files are generated before runninng parallel GRIDSS jobs.",
            ),
            ToolInput(
                "configuration",
                File(optional=True),
                prefix="--configuration",
                doc="configuration file use to override default GRIDSS settings.",
            ),
            ToolInput(
                "labels",
                Array(String, optional=True),
                prefix="--labels",
                separator=",",
                prefix_applies_to_all_elements=False,
                doc='comma separated labels to use in the output VCF for the input files. Supporting read counts for input files with the same label are aggregated (useful for multiple sequencing runs of the same sample). Labels default to input filenames, unless a single read group with a non-empty sample name exists in which case the read group sample name is used (which can be disabled by "useReadGroupSampleNameCategoryLabel=false" in the configuration file). If labels are specified, they must be specified for all input files.',
            ),
            ToolInput(
                "externalaligner",
                String(optional=True),
                prefix="--externalaligner",
                doc="use the system version of bwa instead of the in-process version packaged with GRIDSS",
            ),
            ToolInput(
                "jvmheap",
                String(optional=True),
                prefix="--jvmheap",
                doc="size of JVM heap for assembly and variant calling. (Default: 30g)",
            ),
            ToolInput(
                "maxcoverage",
                Int(optional=True),
                prefix="--maxcoverage",
                doc="maximum coverage. Regions with coverage in excess of this are ignored. (Default: 50000)",
            ),
            ToolInput(
                "picardoptions",
                String(optional=True),
                prefix="--picardoptions",
                doc="additional standard Picard command line options. Useful options include VALIDATION_STRINGENCY=LENIENT and COMPRESSION_LEVEL=0. See https://broadinstitute.github.io/picard/command-line-overview.html",
            ),
            ToolInput(
                "useproperpair",
                String(optional=True),
                prefix="--useproperpair",
                doc="use SAM 'proper pair' flag to determine whether a read pair is discordant. Default: use library fragment size distribution to determine read pair concordance",
            ),
            ToolInput(
                "concordantreadpairdistribution",
                Float(optional=True),
                prefix="--concordantreadpairdistribution",
                doc="portion of 6 sigma read pairs distribution considered concordantly mapped. (Default: 0.995)",
            ),
            ToolInput(
                "keepTempFiles",
                Boolean(optional=True),
                prefix="--keepTempFiles",
                doc="keep intermediate files. Not recommended except for debugging due to the high disk usage.",
            ),
            ToolInput(
                "nojni",
                Boolean(optional=True),
                prefix="--nojni",
                doc="do not use JNI native code acceleration libraries (snappy, GKL, ssw, bwa).",
            ),
            ToolInput(
                "jobindex",
                Int(optional=True),
                prefix="--jobindex",
                doc="zero-based assembly job index (only required when performing parallel assembly across multiple computers)",
            ),
            ToolInput(
                "jobnodes",
                Int(optional=True),
                prefix="--jobnodes",
                doc="total number of assembly jobs (only required when performing parallel assembly across multiple computers). Note than an assembly job with any --job argument is required to be run after all indexed jobs have been completed to gather the output files together.",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput("vcf", Vcf(), glob=InputSelector("outputFilename")),
            ToolOutput(
                "assembly",
                Bam(),
                glob=InputSelector("assemblyFilename"),
            ),
        ]

    def cpus(self, hints: Dict[str, Any]):
        return 8

    def memory(self, hints: Dict[str, Any]):
        return 31

    def bind_metadata(self):

        self.metadata.contributors = ["Jiaan Yu"]
        self.metadata.dateCreated = date(2021, 3, 30)
        self.metadata.dateUpdated = date(2021, 3, 30)
        self.metadata.documentationUrl = (
            "https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation"
        )
        self.metadata.doi = "10.1101/gr.222109.117"
        self.metadata.citation = (
            "Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, "
            "Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss. "
            "GRIDSS: sensitive and specific genomic rearrangement detection using positional "
            "de Bruijn graph assembly. Genome Research, 2017 doi: 10.1101/gr.222109.117"
        )
        self.metadata.documentation = """\
GRIDSS: the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
any step (such as split read identification) can be replaced by another implementation that also outputs
using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.
"""
