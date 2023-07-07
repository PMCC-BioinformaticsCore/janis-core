

from abc import ABC
from typing import List, Any, Dict

from janis_core import (
    ToolOutput,
    ToolInput,
    ToolArgument,
    CommandTool,
    Boolean,
    String,
    File,
    InputSelector,
    CaptureType,
    StringFormatter,
    ToolMetadata,
    CpuSelector,
    MemorySelector,
)
from janis_core import get_value_for_hints_and_ordered_resource_tuple
from ..types import Tsv, FastaWithDict, VcfTabix, BamBai, BedTabix
from .bioinformaticstool import BioinformaticsTool


class IlluminaToolBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "Illumina"

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 16,
            CaptureType.THIRTYX: 32,
            CaptureType.NINETYX: 40,
            CaptureType.THREEHUNDREDX: 40,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 4,
            CaptureType.CHROMOSOME: 32,
            CaptureType.EXOME: 32,
            CaptureType.THIRTYX: 64,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class StrelkaGermlineBase(IlluminaToolBase, ABC):
    def tool(self):
        return "strelka_germline"

    def friendly_name(self):
        return "Strelka (Germline)"

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
        return [
            ToolInput(
                "config",
                File(optional=True),
                prefix="--config",
                position=1,
                shell_quote=False,
                doc="provide a configuration file to override defaults in \
                global config file \
                (/opt/strelka/bin/configureStrelkaGermlineWorkflow.py.ini)",
            ),
            ToolInput(
                "bam",
                BamBai(),
                prefix="--bam",
                position=1,
                shell_quote=False,
                doc="Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated "
                "as each BAM file representing a different sample. [required] (no default)",
            ),
            ToolInput(
                "reference",
                FastaWithDict(),
                prefix="--referenceFasta",
                position=1,
                shell_quote=False,
                doc="samtools-indexed reference fasta file [required]",
            ),
            ToolInput(
                "relativeStrelkaDirectory",
                String(optional=True),
                default="strelka_dir",
                prefix="--runDir",
                position=1,
                shell_quote=False,
                doc="Name of directory to be created where all workflow scripts and output will be written. "
                "Each analysis requires a separate directory.",
            ),
            ToolInput(
                "ploidy",
                VcfTabix(optional=True),
                prefix="--ploidy",
                position=1,
                shell_quote=False,
                doc="Provide ploidy file in VCF. The VCF should include one sample column per input sample "
                "labeled with the same sample names found in the input BAM/CRAM RG header sections. "
                "Ploidy should be provided in records using the FORMAT/CN field, which are interpreted "
                "to span the range [POS+1, INFO/END]. Any CN value besides 1 or 0 will be treated as 2. "
                "File must be tabix indexed. (no default)",
            ),
            ToolInput(
                "noCompress",
                VcfTabix(optional=True),
                prefix="--noCompress",
                position=1,
                shell_quote=False,
                doc="Provide BED file of regions where gVCF block compression is not allowed. "
                "File must be bgzip- compressed/tabix-indexed. (no default)",
            ),
            ToolInput(
                "callContinuousVf",
                String(optional=True),
                prefix="--callContinuousVf",
                doc="Call variants on CHROM without a ploidy prior assumption, "
                "issuing calls with continuous variant frequencies (no default)",
            ),
            ToolInput(
                "rna",
                Boolean(optional=True),
                prefix="--rna",
                position=1,
                shell_quote=False,
                doc="Set options for RNA-Seq input.",
            ),
            ToolInput(
                "indelCandidates",
                VcfTabix(optional=True),
                prefix="--indelCandidates",
                position=1,
                shell_quote=False,
                doc="Specify a VCF of candidate indel alleles. These alleles are always evaluated but only "
                "reported in the output when they are inferred to exist in the sample. "
                "The VCF must be tabix indexed. All indel alleles must be left-shifted/normalized, "
                "any unnormalized alleles will be ignored. This option may be specified more than once, "
                "multiple input VCFs will be merged. (default: None)",
            ),
            ToolInput(
                "forcedGT",
                VcfTabix(optional=True),
                prefix="--forcedGT",
                position=1,
                shell_quote=False,
                doc="Specify a VCF of candidate alleles. These alleles are always evaluated and reported even "
                "if they are unlikely to exist in the sample. The VCF must be tabix indexed. "
                "All indel alleles must be left- shifted/normalized, any unnormalized allele will "
                "trigger a runtime error. This option may be specified more than once, multiple input "
                "VCFs will be merged. Note that for any SNVs provided in the VCF, the SNV site will "
                "be reported (and for gVCF, excluded from block compression), "
                "but the specific SNV alleles are ignored. (default: None)",
            ),
            ToolInput(
                "exome",
                Boolean(optional=True),
                prefix="--exome",
                position=1,
                shell_quote=False,
                doc="Set options for exome note in particular that this flag turns off high-depth filters",
            ),
            ToolInput(
                "targeted",
                Boolean(optional=True),
                prefix="--exome",
                position=1,
                shell_quote=False,
                doc="Set options for other targeted input: "
                "note in particular that this flag turns off high-depth filters",
            ),
            ToolInput(
                tag="callRegions",
                input_type=BedTabix(optional=True),
                prefix="--callRegions=",
                separate_value_from_prefix=False,
                position=1,
                doc="Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of "
                "regions to call. No VCF output will be provided outside of these regions. "
                "The full genome will still be used to estimate statistics from the input "
                "(such as expected depth per chromosome). Only one BED file may be specified. "
                "(default: call the entire genome)",
            ),
            # ToolInput("version", Boolean(optional=True), prefix="--version", position=3, shell_quote=False,
            #           doc="show program's version number and exit"),
            # ToolInput("help", Boolean(optional=True), prefix="--help", position=3, shell_quote=False,
            #           doc="(-h) show this help message and exit"),
            ToolInput(
                "mode",
                String(optional=True),
                default="local",
                prefix="--mode",
                position=3,
                shell_quote=False,
                doc="(-m MODE)  select run mode (local|sge)",
            ),
            ToolInput(
                "queue",
                String(optional=True),
                prefix="--queue",
                position=3,
                shell_quote=False,
                doc="(-q QUEUE) specify scheduler queue name",
            ),
            # ToolInput("dryRun", Boolean(optional=True), prefix="--dryRun", position=3, shell_quote=False,
            #           doc="dryRun (-d,) workflow code without actually running command-tasks"),
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
                "mailTo",
                String(optional=True),
                prefix="--mailTo",
                position=3,
                shell_quote=False,
                doc="(-e) send email notification of job completion status to this address "
                "(may be provided multiple times for more than one email address)",
            ),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput(
                "configPickle",
                File(),
                glob=InputSelector("relativeStrelkaDirectory")
                + "/runWorkflow.py.config.pickle",
            ),
            ToolOutput(
                "script",
                File(),
                glob=InputSelector("relativeStrelkaDirectory") + "/runWorkflow.py",
            ),
            ToolOutput(
                "stats",
                Tsv(),
                glob=InputSelector("relativeStrelkaDirectory")
                + "/results/stats/runStats.tsv",
                doc="A tab-delimited report of various internal statistics from the variant calling process: "
                "Runtime information accumulated for each genome segment, excluding auxiliary steps such "
                "as BAM indexing and vcf merging. Indel candidacy statistics",
            ),
            ToolOutput(
                "variants",
                VcfTabix(),
                glob=InputSelector("relativeStrelkaDirectory")
                + "/results/variants/variants.vcf.gz",
                doc="Primary variant inferences are provided as a series of VCF 4.1 files",
            ),
            ToolOutput(
                "genome",
                VcfTabix(),
                glob=InputSelector("relativeStrelkaDirectory")
                + "/results/variants/genome.vcf.gz",
            ),
        ]

    def arguments(self) -> List[ToolArgument]:
        return [
            ToolArgument(
                "configureStrelkaGermlineWorkflow.py", position=0, shell_quote=False
            ),
            ToolArgument(
                StringFormatter(";")
                + InputSelector("relativeStrelkaDirectory")
                + "/runWorkflow.py",
                position=2,
                shell_quote=False,
            ),
            ToolArgument(
                CpuSelector(None),
                prefix="--jobs",
                position=3,
                shell_quote=False,
                doc=" (-j JOBS)  number of jobs, must be an integer or 'unlimited' "
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

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Illumina",
            doi=None,
            citation=None,  # find citation
            keywords=["broad", "igvtools", "index"],
            documentationUrl="https://github.com/Illumina/strelka",
            documentation="""
Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation 
in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs 
an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively 
selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. 
The germline caller also analyzes input sequencing data using a mixture-model indel error estimation 
method to improve robustness to indel noise. The somatic calling model improves on the original 
Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell 
contamination in the normal sample. A final empirical variant re-scoring step using random forest 
models trained on various call quality features has been added to both callers to further improve precision.

Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average 
indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, 
than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline 
analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis

Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call 
alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant 
reporting uses the gVCF conventions to represent both variant and reference call confidence. 
For best somatic indel performance, Strelka is designed to be run with the Manta structural variant 
and indel caller, which provides additional indel candidates up to a given maxiumum indel size 
(49 by default). By design, Manta and Strelka run together with default settings provide complete 
coverage over all indel sizes (in additional to SVs and SNVs). 

See the user guide for a full description of capabilities and limitations""".strip(),
        )



class Strelka_2_9_9(CommandTool):
    def tool_provider(self):
        return "Illumina"

    def container(self):
        return ""

    def version(self):
        return "2.9.9"


class Strelka_2_9_10(CommandTool, ABC):
    def tool_provider(self):
        return "Illumina"

    def container(self):
        return "michaelfranklin/strelka:2.9.10"

    def version(self):
        return "2.9.10"


class StrelkaGermline_2_9_9(Strelka_2_9_9, StrelkaGermlineBase):
    pass


class StrelkaGermline_2_9_10(Strelka_2_9_10, StrelkaGermlineBase):
    pass


StrelkaGermlineLatest = StrelkaGermline_2_9_10
