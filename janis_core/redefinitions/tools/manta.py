

from abc import ABC, abstractmethod
from typing import Dict, Any, List

from janis_core import (
    Boolean,
    CaptureType,
    CpuSelector,
    File,
    Filename,
    InputSelector,
    Int,
    MemorySelector,
    String,
    StringFormatter,
    ToolArgument,
    ToolInput,
    ToolMetadata,
    ToolOutput,
    get_value_for_hints_and_ordered_resource_tuple,
)

from ..utils import cast_input_bams_to_crams
from ..types import Tsv, BamBai, BedTabix, FastaFai, VcfTabix
from .bioinformaticstool import BioinformaticsTool


class IlluminaToolBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "Illumina"
    

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.CHROMOSOME: 8,
            CaptureType.EXOME: 8,
            CaptureType.THIRTYX: 32,
            CaptureType.NINETYX: 40,
            CaptureType.THREEHUNDREDX: 80,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 8,
            CaptureType.CHROMOSOME: 32,
            CaptureType.EXOME: 32,
            CaptureType.THIRTYX: 64,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class MantaBase(IlluminaToolBase, ABC):
    def tool(self):
        return "manta"

    def base_command(self):
        return None

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 4

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 4

    def inputs(self) -> List[ToolInput]:
        return [*self.config_inputs, *self.running_inputs]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput(
                "python", File(), glob=InputSelector("runDir") + "/runWorkflow.py"
            ),
            ToolOutput(
                "pickle",
                File(),
                glob=InputSelector("runDir") + "/runWorkflow.py.config.pickle",
            ),
            ToolOutput(
                "candidateSV",
                VcfTabix(),
                glob=InputSelector("runDir") + "/results/variants/candidateSV.vcf.gz",
            ),
            ToolOutput(
                "candidateSmallIndels",
                VcfTabix(),
                glob=InputSelector("runDir")
                + "/results/variants/candidateSmallIndels.vcf.gz",
            ),
            ToolOutput(
                "diploidSV",
                VcfTabix(),
                glob=InputSelector("runDir") + "/results/variants/diploidSV.vcf.gz",
            ),
            ToolOutput(
                "alignmentStatsSummary",
                File(),
                glob=InputSelector("runDir")
                + "/results/stats/alignmentStatsSummary.txt",
            ),
            ToolOutput(
                "svCandidateGenerationStats",
                Tsv(),
                glob=InputSelector("runDir")
                + "/results/stats/svCandidateGenerationStats.tsv",
            ),
            ToolOutput(
                "svLocusGraphStats",
                Tsv(),
                glob=InputSelector("runDir") + "/results/stats/svLocusGraphStats.tsv",
            ),
            # optional outputs
            ToolOutput(
                "somaticSV",
                VcfTabix(optional=True),
                glob=InputSelector("runDir") + "/results/variants/somaticSV.vcf.gz",
            ),
            ToolOutput(
                "tumorSV",
                VcfTabix(optional=True),
                glob=InputSelector("runDir") + "/results/variants/tumorSV.vcf.gz",
            ),
        ]

    def arguments(self) -> List[ToolArgument]:
        return [
            ToolArgument("configManta.py", position=0, shell_quote=False),
            ToolArgument(
                StringFormatter(";") + InputSelector("runDir") + "/runWorkflow.py",
                position=2,
                shell_quote=False,
            ),
            ToolArgument(
                CpuSelector(None),
                position=3,
                shell_quote=False,
                prefix="-j",
                doc="(-j) number of jobs, must be an integer or 'unlimited' "
                "(default: Estimate total cores on this node for local mode, 128 for sge mode)",
            ),
            ToolArgument(
                MemorySelector(),
                prefix="--memGb",
                position=3,
                shell_quote=False,
                doc=" (-g MEMGB) gigabytes of memory available to run workflow "
                "-- only meaningful in local mode, must be an integer (default: Estimate the total "
                "memory for this node for local mode, 'unlimited' for sge mode)",
            ),
        ]

    @abstractmethod
    def container(self):
        raise Exception("Strelka version must override docker command")

    def friendly_name(self):
        return "Manta"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2019, 2, 12),
            dateUpdated=date(2019, 2, 19),
            institution="Illumina",
            doi=" doi:10.1093/bioinformatics/btv710",
            citation="Chen, X. et al. (2016) Manta: rapid detection of structural variants and indels for germline and "
            "cancer sequencing applications. Bioinformatics, 32, 1220-1222. doi:10.1093/bioinformatics/btv710",
            keywords=["illumina", "manta", "variant caller"],
            documentationUrl="https://github.com/Illumina/manta",
            documentation="""
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic
variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs,
medium-sized indels and large insertions within a single efficient workflow. The method is
designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is
analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses
can be completed within 2 hours. Manta combines paired and split-read evidence during SV
discovery and scoring to improve accuracy, but does not require split-reads or successful
breakpoint assemblies to report a variant in cases where there is strong evidence otherwise.

It provides scoring models for germline variants in small sets of diploid samples and somatic
variants in matched tumor/normal sample pairs. There is experimental support for analysis of
unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and
reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description
of capabilities and limitations.""".strip(),
        )

    config_inputs = [
        ToolInput(
            "config",
            File(optional=True),
            prefix="--config",
            position=1,
            shell_quote=False,
            doc="provide a configuration file to override defaults in global config file "
            "(/opt/conda/share/manta-1.2.1-0/bin/configManta.py.ini)",
        ),
        ToolInput(
            "bam",
            BamBai(),
            prefix="--bam",
            position=1,
            shell_quote=False,
            doc="FILE Normal sample BAM or CRAM file. May be specified more than once, multiple inputs "
            "will be treated as each BAM file representing a different sample. [optional] (no default)",
        ),
        ToolInput(
            "runDir",
            Filename(),
            prefix="--runDir",
            position=1,
            shell_quote=False,
            doc="Run script and run output will be written to this directory [required] "
            "(default: MantaWorkflow)",
        ),
        ToolInput(
            "reference",
            FastaFai(),
            prefix="--referenceFasta",
            position=1,
            shell_quote=False,
            doc="samtools-indexed reference fasta file [required]",
        ),
        ToolInput(
            "tumorBam",
            BamBai(optional=True),
            prefix="--tumorBam",
            position=1,
            shell_quote=False,
            doc="Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional=null]",
        ),
        ToolInput(
            "exome",
            Boolean(optional=True),
            prefix="--exome",
            position=1,
            shell_quote=False,
            doc="Set options for WES input: turn off depth filters",
        ),
        ToolInput(
            "rna",
            Boolean(optional=True),
            prefix="--rna",
            position=1,
            shell_quote=False,
            doc="Set options for RNA-Seq input. Must specify exactly one bam input file",
        ),
        ToolInput(
            "unstrandedRNA",
            Boolean(optional=True),
            prefix="--unstrandedRNA",
            position=1,
            shell_quote=False,
            doc="Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand",
        ),
        ToolInput(
            "outputContig",
            Boolean(optional=True),
            prefix="--outputContig",
            position=1,
            shell_quote=False,
            doc="Output assembled contig sequences in VCF file",
        ),
        ToolInput(
            "callRegions",
            BedTabix(optional=True),
            prefix="--callRegions",
            position=1,
            shell_quote=False,
            doc="Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of "
            "regions to call. No VCF output will be provided outside of these regions. The full "
            "genome will still be used to estimate statistics from the input (such as expected depth "
            "per chromosome). Only one BED file may be specified. (default: call the entire genome)",
        ),
    ]

    running_inputs = [
        ToolInput(
            "mode",
            String(optional=True),
            default="local",
            prefix="--mode",
            position=3,
            shell_quote=False,
            doc="(-m) select run mode (local|sge)",
        ),
        ToolInput(
            "quiet",
            Boolean(optional=True),
            prefix="--quiet",
            position=3,
            shell_quote=False,
            doc="Don't write any log output to stderr "
            "(but still write to workspace/pyflow.data/logs/pyflow_log.txt)",
        ),
        ToolInput(
            "queue",
            String(optional=True),
            prefix="--queue",
            position=3,
            shell_quote=False,
            doc="(-q) specify scheduler queue name",
        ),
        # ToolInput("dryRun", Boolean(optional=True), prefix="--dryRun", position=3, shell_quote=False,
        #           doc="(-d) dryRun workflow code without actually running command - tasks"),
        ToolInput(
            "maxTaskRuntime",
            String(optional=True),
            prefix="--maxTaskRuntime",
            position=3,
            shell_quote=False,
            doc="(format: hh:mm:ss) Specify scheduler max runtime per task, argument is "
            "provided to the 'h_rt' resource limit if using SGE (no default)",
        ),
    ]



class MantaCramBase(MantaBase):
    def id(self):
        return super().id() + "_cram"

    def inputs(self):
        # we want every input which is a bam in the original to be a cram now
        return cast_input_bams_to_crams(super().inputs())



class Manta_1_4_0(MantaBase):
    def container(self):
        return "michaelfranklin/manta:1.4.0"  # add 1.4.0 after it's released

    def version(self):
        return "1.4.0"


class Manta_1_5_0(MantaBase):
    def container(self):
        return "michaelfranklin/manta:1.5.0"  # add 1.4.0 after it's released

    def version(self):
        return "1.5.0"


class MantaCram_1_4_0(MantaCramBase):
    def container(self):
        return "michaelfranklin/manta:1.4.0"  # add 1.4.0 after it's released

    def version(self):
        return "1.4.0"


class MantaCram_1_5_0(MantaCramBase):
    def container(self):
        return "michaelfranklin/manta:1.5.0"  # add 1.4.0 after it's released

    def version(self):
        return "1.5.0"


MantaLatest = Manta_1_5_0
