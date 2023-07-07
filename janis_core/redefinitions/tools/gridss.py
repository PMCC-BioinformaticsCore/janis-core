
from datetime import date
from typing import List, Dict, Any

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
    CaptureType,
    CpuSelector,
    get_value_for_hints_and_ordered_resource_tuple,
)
from ..types import Bam, BamBai, FastaWithDict, Bed, Vcf
from .bioinformaticstool import BioinformaticsTool


### GRIDSS 2.2 ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 8,
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
            # https://github.com/PapenfussLab/gridss#how-much-memory-should-i-give-gridss
            CaptureType.TARGETED: 31,
            CaptureType.CHROMOSOME: 31,
            CaptureType.EXOME: 31,
            CaptureType.THIRTYX: 31,
            CaptureType.NINETYX: 31,
            CaptureType.THREEHUNDREDX: 31,
        },
    )
]

class GridssBase_2_2(BioinformaticsTool):
    
    def tool(self) -> str:
        return "gridss"

    def tool_provider(self):
        return "Papenfuss Labs"

    def friendly_name(self) -> str:
        return "Gridss"

    def base_command(self):
        return [
            "java",
            "-XX:+UnlockExperimentalVMOptions",
            "-XX:+UseCGroupMemoryLimitForHeap",
            "-XX:MaxRAMFraction=1",
            "-XshowSettings:vm",
            "-Dsamjdk.create_index=true",
            "-Dsamjdk.use_async_io_read_samtools=true",
            "-Dsamjdk.use_async_io_write_samtools=true",
            "-Dsamjdk.use_async_io_write_tribble=true",
            "-Dgridss.gridss.output_to_temp_file=true",
            "-Dsamjdk.buffer_size=4194304",
            "-cp",
            # insert the following line within the VersionedTool:
            # /data/gridss/gridss-$version-gridss-jar-with-dependencies.jar gridss.CallVariants
        ]

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf"),
                prefix="OUTPUT=",
                separate_value_from_prefix=False,
                doc="(O=) VCF structural variation calls. Required.",
            ),
            ToolInput(
                "reference",
                FastaWithDict(),
                prefix="REFERENCE_SEQUENCE=",
                separate_value_from_prefix=False,
            ),
            ToolInput(
                "bams",
                Array(BamBai()),
                prefix="INPUT=",
                separate_value_from_prefix=False,
                prefix_applies_to_all_elements=True,
                doc="(I=File Coordinate-sorted input BAM file. Default value: null. "
                "This option may be specified 0 or more times.",
            ),
            ToolInput(
                "assemblyFilename",
                Filename(suffix=".assembled", extension=".bam"),
                prefix="ASSEMBLY=",
                separate_value_from_prefix=False,
                doc="Breakend assemblies which have undergone split read identification Required.",
            ),
            ToolInput(
                "inputLabel",
                String(optional=True),
                prefix="INPUT_LABEL=",
                separate_value_from_prefix=False,
                doc="Input label. Variant calling evidence breakdowns are reported for each label. Default "
                "labels correspond to INPUT filenames. When specifying labels, labels must be provided for "
                "all input files. Default value: null. This option may be specified 0 or more times.",
            ),
            ToolInput(
                "inputMaxFragmentSize",
                Int(optional=True),
                prefix="INPUT_MAX_FRAGMENT_SIZE=",
                separate_value_from_prefix=False,
                doc="Per input maximum concordant fragment size. Default value: null. "
                "This option may be specified 0 or more times.",
            ),
            ToolInput(
                "inputMinFragmentSize",
                Int(optional=True),
                prefix="INPUT_MIN_FRAGMENT_SIZE=",
                separate_value_from_prefix=False,
                doc="Per input minimum concordant fragment size. Default value: null. "
                "This option may be specified 0 or more times.",
            ),
            ToolInput(
                "readPairConcordantPercent",
                Float(optional=True),
                prefix="READ_PAIR_CONCORDANT_PERCENT=",
                separate_value_from_prefix=False,
                doc="Percent of read pairs considered concorant (0.0-1.0). If this is unset, the SAM proper "
                "pair flag is used to determine whether a read is discordantly aligned. Explicit fragment "
                "size specification overrides this setting. Default value: 0.995. "
                "This option can be set to 'null' to clear the default value.",
            ),
            ToolInput(
                "blacklist",
                Bed(optional=True),
                prefix="BLACKLIST=",
                separate_value_from_prefix=False,
                doc="(BL=File) BED blacklist of regions to ignore. Assembly of regions such as high-coverage "
                "centromeric repeats is slow, and if such regions are to be filtered in downstream "
                "analysis anyway, blacklisting those region will improve runtime performance. "
                "For human WGS, the ENCODE DAC blacklist is recommended. Default value: null.",
            ),
            ToolInput(
                "configurationFile",
                File(optional=True),
                prefix="CONFIGURATION_FILE=",
                separate_value_from_prefix=False,
                doc="(C=File) gridss configuration file containing overrides Default value: null.",
            ),
            ToolInput(
                "workerThreads",
                Int(optional=True),
                prefix="WORKER_THREADS=",
                separate_value_from_prefix=False,
                doc="(THREADS=Integer  Number of worker threads to spawn. Defaults to number of cores available. "
                "Note that I/O threads are not included in this worker thread count so CPU usage can be "
                "higher than the number of worker thread. Default value: 6. "
                "This option can be set to 'null' to clear the default value.",
            ),
            ToolInput(
                "workingDir",
                String(optional=True),
                prefix="WORKING_DIR=",
                default=".",
                separate_value_from_prefix=False,
                doc="Directory to place intermediate results directories. Default location is the same "
                "directory as the associated input or output file. Default value: null.",
            ),
            ToolInput(
                "ignoreDuplicates",
                Boolean(optional=True),
                prefix="IGNORE_DUPLICATES=",
                separate_value_from_prefix=False,
                doc="Ignore reads marked as duplicates. Default value: true. This option can be set to 'null' "
                "to clear the default value. Possible values: {true, false}",
            ),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("vcf", Vcf(), glob=InputSelector("outputFilename")),
            ToolOutput(
                "assembly",
                BamBai(),
                glob=InputSelector("assemblyFilename"),
                secondaries_present_as={".bai": "^.bai"},
            ),
        ]

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 8

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 31

    def bind_metadata(self):

        self.metadata.contributors = ["Michael Franklin"]
        self.metadata.dateCreated = date(2019, 6, 19)
        self.metadata.dateUpdated = date(2019, 7, 3)
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


### GRIDSS 2.4 ###


CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 8,
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
            # https://github.com/PapenfussLab/gridss#how-much-memory-should-i-give-gridss
            CaptureType.TARGETED: 31,
            CaptureType.CHROMOSOME: 31,
            CaptureType.EXOME: 31,
            CaptureType.THIRTYX: 31,
            CaptureType.NINETYX: 31,
            CaptureType.THREEHUNDREDX: 31,
        },
    )
]

class GridssBase_2_4(BioinformaticsTool):
    def tool(self) -> str:
        return "gridss"

    def tool_provider(self):
        return "Papenfuss Labs"

    def friendly_name(self) -> str:
        return "Gridss"

    def base_command(self):
        return "/opt/gridss/gridss.sh"

    def inputs(self):
        return [
            ToolInput("bams", Array(BamBai()), position=10),
            ToolInput("reference", FastaWithDict(), position=1, prefix="--reference"),
            ToolInput(
                "outputFilename",
                Filename(suffix=".svs", extension=".vcf"),
                position=2,
                prefix="--output",
            ),
            ToolInput(
                "assemblyFilename",
                Filename(suffix=".assembled", extension=".bam"),
                position=3,
                prefix="--assembly",
            ),
            ToolInput(
                "threads", Int(optional=True), default=CpuSelector(), prefix="--threads"
            ),
            ToolInput(
                "blacklist", Bed(optional=True), position=4, prefix="--blacklist"
            ),
            ToolInput(
                "tmpdir", String(optional=True), default="./TMP", prefix="--workingdir"
            ),
        ]

    def outputs(self):
        return [
            ToolOutput("out", Vcf(), glob=InputSelector("outputFilename")),
            ToolOutput(
                "assembly",
                Bam(),
                glob=InputSelector("assemblyFilename"),
            ),
        ]

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 8

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 31

    def bind_metadata(self):

        self.metadata.contributors = ["Michael Franklin"]
        self.metadata.dateCreated = date(2019, 6, 19)
        self.metadata.dateUpdated = date(2019, 8, 20)
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



### GRIDSS 2.10 ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 8,
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
            # https://github.com/PapenfussLab/gridss#how-much-memory-should-i-give-gridss
            CaptureType.TARGETED: 31,
            CaptureType.CHROMOSOME: 31,
            CaptureType.EXOME: 31,
            CaptureType.THIRTYX: 31,
            CaptureType.NINETYX: 31,
            CaptureType.THREEHUNDREDX: 31,
        },
    )
]

class GridssBase_2_10(BioinformaticsTool):
    def tool(self) -> str:
        return "gridss"

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
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 8

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
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



### OTHER VERSIONS ###

class Gridss_2_2_3(GridssBase_2_2):
    def base_command(self):
        return [
            *super().base_command(),
            "/data/gridss/gridss-2.2.3-gridss-jar-with-dependencies.jar",
            "gridss.CallVariants",
        ]

    def container(self):
        return "gridss/gridss:v2.2.3"

    def version(self):
        return "v2.2.3"


class Gridss_2_4_0(GridssBase_2_2):
    def base_command(self):
        return [
            *super().base_command(),
            "/data/gridss/gridss-2.4.0-gridss-jar-with-dependencies.jar",
            "gridss.CallVariants",
        ]

    def container(self):
        return "gridss/gridss:2.4.0"

    def version(self):
        return "v2.4.0"


class Gridss_2_5_1(GridssBase_2_4):
    def base_command(self):
        return "gridss.sh"

    def container(self):
        return "michaelfranklin/gridss:2.5.1-dev2"

    def version(self):
        return "v2.5.1-dev"


class Gridss_2_6_2(GridssBase_2_4):
    def container(self):
        # https://hub.docker.com/r/gridss/gridss
        return "gridss/gridss:2.6.2"

    def version(self):
        return "v2.6.2"


class Gridss_2_9_4(GridssBase_2_4):
    def container(self):
        return "gridss/gridss:2.9.4"

    def version(self) -> str:
        return "v2.9.4"


# 2.8.3 is last version before library optimisation
class Gridss_2_8_3(GridssBase_2_10):
    def container(self):
        return "gridss/gridss:2.8.3"

    def version(self) -> str:
        return "v2.8.3"


class Gridss_2_10_2(GridssBase_2_10):
    def container(self):
        return "gridss/gridss:2.10.2"

    def version(self) -> str:
        return "v2.10.2"


GridssLatest = Gridss_2_10_2
