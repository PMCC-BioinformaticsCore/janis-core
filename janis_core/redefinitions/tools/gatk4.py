
import os
import operator
from datetime import datetime

from abc import ABC, abstractmethod
from typing import Dict, Any, Union, List

from janis_core import (
    ToolInput,
    Filename,
    String,
    Array,
    File,
    Int,
    Boolean,
    ToolOutput,
    InputSelector,
    ToolArgument,
    CaptureType,
    StringFormatter,
    MemorySelector,
    ToolMetadata,
    Double,
    Float
)
from janis_core.operators.logical import If, IsDefined
from janis_core.operators.standard import JoinOperator, FirstOperator
from janis_core import get_value_for_hints_and_ordered_resource_tuple

from ..types import Bam, BamBai, FastaWithDict, TextFile, Tsv, Bed, VcfTabix
from ..tools import BioinformaticsTool

from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)



### BASE ###

class Gatk4ToolBase(BioinformaticsTool, ABC):

    def tool_provider(self):
        return "GATK4"

    @classmethod
    def base_command(cls):
        return ["gatk", cls.gatk_command()]

    @classmethod
    @abstractmethod
    def gatk_command(cls):
        raise Exception("Subclass must override 'gatk_command' method")

    def inputs(self):
        return [
            ToolInput("javaOptions", Array(String, optional=True)),
            ToolInput(
                "compression_level",
                Int(optional=True),
                doc="Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.",
            )
            # ToolInput("pg-tag", Boolean(optional=True), prefix="--add-output-sam-program-record",
            #           doc="If true, adds a PG tag to created SAM/BAM/CRAM files.")
        ]

    @abstractmethod
    def container(self):
        raise Exception(
            "An error likely occurred when resolving the method order for docker for the Gatk classes "
            "or you're trying to execute the docker method of the base class (ie, don't do that). "
            "The method order resolution must preference Gatkbase subclasses, "
            "and the subclass must contain a definition for docker."
        )

    def arguments(self):
        return [
            ToolArgument(
                StringFormatter(
                    "-Xmx{memory}G {compression} {otherargs}",
                    memory=MemorySelector() * 3 / 4,
                    compression=If(
                        IsDefined(InputSelector("compression_level")),
                        "-Dsamjdk.compress_level=" + InputSelector("compression_level"),
                        "",
                    ),
                    otherargs=JoinOperator(
                        FirstOperator([InputSelector("javaOptions"), []]), " "
                    ),
                ),
                prefix="--java-options",
                position=-1,
            )
        ]



### VERSIONS ###

class Gatk_4_0_12(ABC):
    def container(self):
        return "broadinstitute/gatk:4.0.12.0"

    def version(self):
        return "4.0.12.0"


class Gatk_4_1_2_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.2.0"

    def version(self):
        return "4.1.2.0"


class Gatk_4_1_3_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.3.0"

    def version(self):
        return "4.1.3.0"


class Gatk_4_1_4_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.4.0"

    def version(self):
        return "4.1.4.0"


class Gatk_4_1_4_1(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.4.1"

    def version(self):
        return "4.1.4.1"


class Gatk_4_1_5_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.5.0"

    def version(self):
        return "4.1.5.0"


class Gatk_4_1_6_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.6.0"

    def version(self):
        return "4.1.6.0"


class Gatk_4_1_7_0(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.7.0"

    def version(self):
        return "4.1.7.0"


class Gatk_4_1_8_1(ABC):
    def container(self):
        return "broadinstitute/gatk:4.1.8.1"

    def version(self):
        return "4.1.8.1"


Gatk4Latest = Gatk_4_1_8_1



### HAPLOTYPE CALLER ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 32,
            CaptureType.EXOME: 32,
            CaptureType.THIRTYX: 32,
            CaptureType.NINETYX: 32,
            CaptureType.THREEHUNDREDX: 32,
        },
    )
]


class Gatk4HaplotypeCallerBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "HaplotypeCaller"

    def tool(self):
        return "Gatk4HaplotypeCaller"

    def friendly_name(self):
        return "GATK4: Haplotype Caller"

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

    def inputs(self):
        return [
            *super(Gatk4HaplotypeCallerBase, self).inputs(),
            *Gatk4HaplotypeCallerBase.optional_args,
            ToolInput(
                "inputRead",
                BamBai(),
                doc="BAM/SAM/CRAM file containing reads",
                prefix="--input",
                secondaries_present_as={".bai": "^.bai"},
            ),
            ToolInput(
                "reference",
                FastaWithDict(),
                position=5,
                prefix="--reference",
                doc="Reference sequence file",
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("inputRead", remove_file_extension=True),
                    extension=".vcf.gz",
                ),
                position=8,
                prefix="--output",
                doc="File to which variants should be written",
            ),
            ToolInput(
                "dbsnp",
                VcfTabix(optional=True),
                position=7,
                prefix="--dbsnp",
                doc="(Also: -D) A dbSNP VCF file.",
            ),
            ToolInput(
                "intervals",
                Bed(optional=True),
                prefix="--intervals",
                doc="-L (BASE) One or more genomic intervals over which to operate",
            ),
            ToolInput(
                "outputBamName",
                Filename(
                    prefix=InputSelector("inputRead", remove_file_extension=True),
                    extension=".bam",
                ),
                position=8,
                prefix="-bamout",
                doc="File to which assembled haplotypes should be written",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                VcfTabix,
                glob=InputSelector("outputFilename"),
                doc="A raw, unfiltered, highly sensitive callset in VCF format. "
                "File to which variants should be written",
            ),
            ToolOutput(
                "bam",
                BamBai,
                glob=InputSelector("outputBamName"),
                doc="File to which assembled haplotypes should be written",
                secondaries_present_as={".bai": "^.bai"},
            ),
        ]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "haplotype"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#",
            documentation="""Call germline SNPs and indels via local re-assembly of haplotypes
    
The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes 
in an active region. In other words, whenever the program encounters a region showing signs of variation, it 
discards the existing mapping information and completely reassembles the reads in that region. This allows the 
HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when 
they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at 
calling indels than position-based callers like UnifiedGenotyper.

In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to 
generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint 
genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of 
samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).

In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele 
frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. 
For that purpose, use Mutect2 instead.

Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge 
for most variant callers, on the condition that the input read data has previously been processed according 
to our recommendations as documented (https://software.broadinstitute.org/gatk/documentation/article?id=4067).
""".strip(),
        )

    optional_args = [
        ToolInput(
            "pairHmmImplementation",
            String(optional=True),
            prefix="--pair-hmm-implementation",
            doc="The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE",
        ),
        ToolInput(
            "activityProfileOut",
            String(optional=True),
            prefix="--activity-profile-out",
            doc="Output the raw activity profile results in IGV format (default: null)",
        ),
        ToolInput(
            "alleles",
            File(optional=True),
            prefix="--alleles",
            doc="(default: null) The set of alleles at which to genotype when --genotyping_mode "
            "is GENOTYPE_GIVEN_ALLELES",
        ),
        ToolInput(
            "annotateWithNumDiscoveredAlleles",
            Boolean(optional=True),
            prefix="--annotate-with-num-discovered-alleles",
            doc="If provided, we will annotate records with the number of alternate alleles that were "
            "discovered (but not necessarily genotyped) at a given site",
        ),
        ToolInput(
            "annotation",
            Array(String(), optional=True),
            prefix="--annotation",
            doc="-A: One or more specific annotations to add to variant calls",
        ),
        ToolInput(
            "annotationGroup",
            Array(String(), optional=True),
            prefix="--annotation-group",
            doc="-G	One or more groups of annotations to apply to variant calls",
        ),
        ToolInput(
            "annotationsToExclude",
            Array(String(), optional=True),
            prefix="--annotations-to-exclude",
            doc="-AX	One or more specific annotations to exclude from variant calls",
        ),
        ToolInput(
            "arguments_file",
            Array(File(), optional=True),
            prefix="--arguments_file",
            doc="read one or more arguments files and add them to the command line",
        ),
        ToolInput(
            "assemblyRegionOut",
            String(optional=True),
            prefix="--assembly-region-out",
            doc="(default: null) Output the assembly region to this IGV formatted file. Which annotations to "
            "exclude from output in the variant calls. Note that this argument has higher priority than "
            "the -A or -G arguments, so these annotations will be excluded even if they are explicitly "
            "included with the other options.",
        ),
        ToolInput(
            "baseQualityScoreThreshold",
            Int(optional=True),
            prefix="--base-quality-score-threshold",
            doc="(default: 18) Base qualities below this threshold will be reduced to the minimum (6)",
        ),
        ToolInput(
            "cloudIndexPrefetchBuffer",
            Int(optional=True),
            prefix="--cloud-index-prefetch-buffer",
            doc="-CIPB (default: -1) Size of the cloud-only prefetch buffer (in MB; 0 to disable). "
            "Defaults to cloudPrefetchBuffer if unset.",
        ),
        ToolInput(
            "cloudPrefetchBuffer",
            Int(optional=True),
            prefix="--cloud-prefetch-buffer",
            doc="-CPB (default: 40) Size of the cloud-only prefetch buffer (in MB; 0 to disable).",
        ),
        ToolInput(
            "contaminationFractionToFilter",
            Double(optional=True),
            prefix="--contamination-fraction-to-filter",
            doc="-contamination (default: 0.0) Fraction of contamination in sequencing data "
            "(for all samples) to aggressively remove",
        ),
        ToolInput(
            "correctOverlappingQuality",
            Boolean(optional=True),
            prefix="--correct-overlapping-quality",
            doc="Undocumented option",
        ),
        # ToolInput("dbsnp", VcfIdx(optional=True), prefix="--dbsnp", doc="-D (default: null) dbSNP file"),
        ToolInput(
            "disableBamIndexCaching",
            Boolean(optional=True),
            prefix="--disable-bam-index-caching",
            doc="-DBIC. If true, don't cache bam indexes, this will reduce memory requirements but may harm "
            "performance if many intervals are specified. Caching is automatically disabled if "
            "there are no intervals specified.",
        ),
        # ToolInput("disableSequenceDictionaryValidation", Boolean(optional=True), prefix="--disable-sequence-dictionary-validation",
        #           doc="If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!"),
        ToolInput(
            "founderId",
            Array(String(), optional=True),
            prefix="--founder-id",
            doc='Samples representing the population "founders"',
        ),
        # ToolInput("gcsMaxRetries", Int(optional=True), prefix="--gcs-max-retries",
        #           doc="-gcs-retries (default: 20) If the GCS bucket channel errors out, "
        #               "how many times it will attempt to re-initiate the connection"),
        # ToolInput("gcsProjectForRequesterPays", String(), prefix="--gcs-project-for-requester-pays",
        #           doc="Project to bill when accessing \"requester pays\" buckets. If unset, these buckets cannot be accessed."),
        ToolInput(
            "genotypingMode",
            String(optional=True),
            prefix="--genotyping-mode",
            doc="(default: DISCOVERY) Specifies how to determine the alternate alleles to use for genotyping. "
            "The --genotyping-mode argument is an enumerated type (GenotypingOutputMode), which can have one "
            "of the following values: DISCOVERY (The genotyper will choose the most likely alternate allele) "
            "or GENOTYPE_GIVEN_ALLELES (Only the alleles passed by the user should be considered).",
        ),
        # ToolInput("graphOutput", DataType(optional=True), prefix="--graph-output", doc="-graph	null	Write debug assembly graph information to this file"),
        ToolInput(
            "heterozygosity",
            Double(optional=True),
            prefix="--heterozygosity",
            doc="(default: 0.001) Heterozygosity value used to compute prior likelihoods for any locus. The "
            "expected heterozygosity value used to compute prior probability that a locus is non-reference. "
            "The default priors are for provided for humans: het = 1e-3 which means that the probability "
            "of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i) Note that heterozygosity as "
            "used here is the population genetics concept: "
            "http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics . "
            "That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population "
            "of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. "
            "Note that this quantity has nothing to do with the likelihood of any given sample having a "
            "heterozygous genotype, which in the GATK is purely determined by the probability of the observed "
            "data P(D | AB) under the model that there may be a AB het genotype. The posterior probability "
            "of this AB genotype would use the het prior, but the GATK only uses this posterior probability "
            "in determining the prob. that a site is polymorphic. So changing the het parameters only "
            "increases the chance that a site will be called non-reference across all samples, but doesn't "
            "actually change the output genotype likelihoods at all, as these aren't posterior probabilities "
            "at all. The quantity that changes whether the GATK considers the possibility of a het genotype "
            "at all is the ploidy, which determines how many chromosomes each individual in the species carries.",
        ),
        ToolInput(
            "heterozygosityStdev",
            Double(optional=True),
            prefix="--heterozygosity-stdev",
            doc="(default 0.01) Standard deviation of heterozygosity for SNP and indel calling.",
        ),
        ToolInput(
            "indelHeterozygosity",
            Double(optional=True),
            prefix="--indel-heterozygosity",
            doc="(default: 1.25E-4) Heterozygosity for indel calling. This argument informs the prior "
            "probability of having an indel at a site. (See heterozygosity)",
        ),
        ToolInput(
            "intervalMergingRule",
            String(optional=True),
            prefix="--interval-merging-rule",
            doc="-imr (default: ALL) Interval merging rule for abutting intervals. By default, the program "
            "merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually "
            "overlap) into a single continuous interval. However you can change this behavior if you want "
            "them to be treated as separate intervals instead. The --interval-merging-rule argument is an "
            "enumerated type (IntervalMergingRule), which can have one of the following values:"
            "[ALL, OVERLAPPING]",
        ),
        ToolInput(
            "maxReadsPerAlignmentStart",
            Int(optional=True),
            prefix="--max-reads-per-alignment-start",
            doc="(default: 50) Maximum number of reads to retain per alignment start position. "
            "Reads above this threshold will be downsampled. Set to 0 to disable.",
        ),
        ToolInput(
            "minBaseQualityScore",
            Int(optional=True),
            prefix="--min-base-quality-score",
            doc="-mbq (default: 10) Minimum base quality required to consider a base for calling",
        ),
        ToolInput(
            "nativePairHmmThreads",
            Int(optional=True),
            prefix="--native-pair-hmm-threads",
            doc="(default: 4) How many threads should a native pairHMM implementation use",
        ),
        ToolInput(
            "nativePairHmmUseDoublePrecision",
            Boolean(optional=True),
            prefix="--native-pair-hmm-use-double-precision",
            doc="use double precision in the native pairHmm. "
            "This is slower but matches the java implementation better",
        ),
        ToolInput(
            "numReferenceSamplesIfNoCall",
            Int(optional=True),
            prefix="--num-reference-samples-if-no-call",
            doc="(default: 0) Number of hom-ref genotypes to infer at sites not present in a panel. When a "
            "variant is not seen in any panel, this argument controls whether to infer (and with what "
            "effective strength) that only reference alleles were observed at that site. "
            'E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000".',
        ),
        ToolInput(
            "outputMode",
            String(optional=True),
            prefix="--output-mode",
            doc="(default: EMIT_VARIANTS_ONLY) Specifies which type of calls we should output. The --output-mode "
            "argument is an enumerated type (OutputMode), which can have one of the following values: "
            "[EMIT_VARIANTS_ONLY (produces calls only at variant sites), "
            "EMIT_ALL_CONFIDENT_SITES (produces calls at variant sites and confident reference sites), "
            "EMIT_ALL_SITES (produces calls at any callable site regardless of confidence; "
            "this argument is intended only for point mutations (SNPs) in DISCOVERY mode or "
            "generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce "
            "a comprehensive set of indels in DISCOVERY mode)]",
        ),
        ToolInput(
            "pedigree",
            File(optional=True),
            prefix="--pedigree",
            doc='-ped (default: null) Pedigree file for determining the population "founders"',
        ),
        ToolInput(
            "populationCallset",
            File(optional=True),
            prefix="--population-callset",
            doc="-population (default: null) Callset to use in calculating genotype priors",
        ),
        ToolInput(
            "sampleName",
            String(optional=True),
            prefix="--sample-name",
            doc="-ALIAS (default: null) Name of single sample to use from a multi-sample bam. You can use this "
            "argument to specify that HC should process a single sample out of a multisample BAM file. "
            "This is especially useful if your samples are all in the same file but you need to run them "
            "individually through HC in -ERC GVC mode (which is the recommended usage). "
            "Note that the name is case-sensitive.",
        ),
        ToolInput(
            "samplePloidy",
            Int(optional=True),
            prefix="--sample-ploidy",
            doc="-ploidy (default: 2) Ploidy (number of chromosomes) per sample. "
            "For pooled data, set to (Number of samples in each pool * Sample Ploidy). "
            "Sample ploidy - equivalent to number of chromosomes per pool. In pooled "
            "experiments this should be = # of samples in pool * individual sample ploidy",
        ),
        ToolInput(
            "sitesOnlyVcfOutput",
            Boolean(optional=True),
            prefix="--sites-only-vcf-output",
            doc="(default: false) If true, don't emit genotype fields when writing vcf file output.",
        ),
        ToolInput(
            "standardMinConfidenceThresholdForCalling",
            Double(optional=True),
            prefix="--standard-min-confidence-threshold-for-calling",
            doc="-stand-call-conf (default: 10.0) The minimum phred-scaled confidence "
            "threshold at which variants should be called",
        ),
        ToolInput(
            "useNewQualCalculator",
            Boolean(optional=True),
            prefix="--use-new-qual-calculator",
            doc="-new-qual If provided, we will use the new AF model instead of the so-called exact model",
        ),
        ToolInput(
            "gvcfGqBands",
            Array(Int, optional=True),
            prefix="-GQB",
            prefix_applies_to_all_elements=True,
            doc="(--gvcf-gq-bands) Exclusive upper bounds for reference confidence GQ"
            " bands (must be in [1, 100] and specified in increasing order)",
        ),
        ToolInput(
            "emitRefConfidence",
            String(optional=True),
            prefix="--emit-ref-confidence",
            doc="(-ERC) Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)",
        ),
        ToolInput(
            "dontUseSoftClippedBases",
            Boolean(optional=True),
            prefix="--dont-use-soft-clipped-bases",
            doc="Do not analyze soft clipped bases in the reads",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "inputRead": f"{remote_dir}/NA12878-BRCA1.split.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                    "dbsnp": f"{remote_dir}/Homo_sapiens_assembly38.dbsnp138.BRCA1.vcf.gz",
                    "javaOptions": ["-Xmx6G"],
                    "pairHmmImplementation": "LOGLESS_CACHING",
                },
                output=VcfTabix.basic_test(
                    "out",
                    12800,
                    270,
                    214,
                    ["GATKCommandLine"],
                    "0224e24e5fc27286ee90c8d3c63373a7",
                )
                + BamBai.basic_test(
                    "bam",
                    596698,
                    21272,
                    f"{remote_dir}/NA12878-BRCA1.haplotyped.flagstat",
                    "d83b4c0d8eab24a3be1cc6af4f827753",
                    "b4bb4028b8679a3a635e3ad87126a097",
                ),
            )
        ]


class Gatk4HaplotypeCaller_4_1_3(Gatk_4_1_3_0, Gatk4HaplotypeCallerBase):
    pass





### Gatk4SplitReads ###

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 8,
            CaptureType.EXOME: 8,
            CaptureType.THIRTYX: 8,
            CaptureType.NINETYX: 16,
            CaptureType.THREEHUNDREDX: 16,
        },
    )
]


class Gatk4SplitReadsBase(Gatk4ToolBase):
    def friendly_name(self) -> str:
        return "GATK4: SplitReads"

    def tool(self) -> str:
        return "Gatk4SplitReads"

    @classmethod
    def gatk_command(cls):
        return "SplitReads"

    def directories_to_create(self) -> Union[str, List[str]]:
        return [InputSelector("outputFilename")]

    def inputs(self):
        return [
            ToolInput(
                "outputFilename",
                String(),
                prefix="--output",
                default=".",
                doc="The directory to output SAM/BAM/CRAM files. Default value: '.' ",
            ),
            ToolInput(
                "bam",
                BamBai,
                prefix="--input",
                position=1,
                secondaries_present_as={".bai": "^.bai"},
                doc="(-I:String) BAM/SAM/CRAM file containing reads  This argument must be specified at least once.",
            ),
            ToolInput(
                tag="intervals",
                input_type=Bed(optional=True),
                prefix="--intervals",
                doc="(-L:String) One or more genomic intervals over which to operate This argument may be specified 0 or more times. Default value: null. ",
            ),
            *super().inputs(),
            *Gatk4SplitReadsBase.additional_args,
        ]

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 4

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai,
                glob=(InputSelector("outputFilename") + '/' + InputSelector("bam").basename()),
                doc="Bam",
                secondaries_present_as={".bai": "^.bai"},
            )
        ]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 9, 16),
            dateUpdated=datetime(2019, 9, 16),
            documentation="USAGE: SplitReads [arguments]\nOutputs reads from a SAM/BAM/CRAM by read group, sample and library name\nVersion:4.1.3.0",
        )

    additional_args = [
        ToolInput(
            tag="addOutputSamProgramRecord",
            input_type=Boolean(optional=True),
            prefix="-add-output-sam-program-record",
            doc="(--add-output-sam-program-record)  If true, adds a PG tag to created SAM/BAM/CRAM files.  Default value: true. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="addOutputVcfCommandLine",
            input_type=Boolean(optional=True),
            prefix="-add-output-vcf-command-line",
            doc="(--add-output-vcf-command-line)  If true, adds a command line header line to created VCF files.  Default value: true. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="arguments_file",
            input_type=File(optional=True),
            prefix="--arguments_file:File",
            doc="read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null. ",
        ),
        ToolInput(
            tag="cloudIndexPrefetchBuffer",
            input_type=String(optional=True),
            prefix="--cloud-index-prefetch-buffer",
            doc="(-CIPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.  Default value: -1. ",
        ),
        ToolInput(
            tag="cloudPrefetchBuffer",
            input_type=String(optional=True),
            prefix="--cloud-prefetch-buffer",
            doc="(-CPB:Integer)  Size of the cloud-only prefetch buffer (in MB; 0 to disable).  Default value: 40. ",
        ),
        ToolInput(
            tag="createOutputBamIndex",
            input_type=String(optional=True),
            prefix="--create-output-bam-index",
            doc="(-OBI:Boolean)  If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM file.  Default value: true. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="createOutputBamMd5",
            input_type=String(optional=True),
            prefix="--create-output-bam-md5",
            doc="(-OBM:Boolean)  If true, create a MD5 digest for any BAM/SAM/CRAM file created  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="createOutputVariantIndex",
            input_type=String(optional=True),
            prefix="--create-output-variant-index",
            doc="(-OVI:Boolean)  If true, create a VCF index when writing a coordinate-sorted VCF file.  Default value: true. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="createOutputVariantMd5",
            input_type=String(optional=True),
            prefix="--create-output-variant-md5",
            doc="(-OVM:Boolean)  If true, create a a MD5 digest any VCF file created.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="disableBamIndexCaching",
            input_type=String(optional=True),
            prefix="--disable-bam-index-caching",
            doc="(-DBIC:Boolean)  If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified.  Caching is automatically disabled if there are no intervals specified.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="disableReadFilter",
            input_type=String(optional=True),
            prefix="--disable-read-filter",
            doc="(-DF:String)  Read filters to be disabled before analysis  This argument may be specified 0 or more times. Default value: null. Possible Values: {WellformedReadFilter}",
        ),
        ToolInput(
            tag="disableSequenceDictionaryValidation",
            input_type=Boolean(optional=True),
            prefix="-disable-sequence-dictionary-validation",
            doc="(--disable-sequence-dictionary-validation)  If specified, do not check the sequence dictionaries from our inputs for compatibility. Use at your own risk!  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="excludeIntervals",
            input_type=String(optional=True),
            prefix="--exclude-intervals",
            doc="(-XL:StringOne) This argument may be specified 0 or more times. Default value: null. ",
        ),
        ToolInput(
            tag="gatkConfigFile",
            input_type=File(optional=True),
            prefix="--gatk-config-file",
            doc="A configuration file to use with the GATK. Default value: null.",
        ),
        ToolInput(
            tag="gcsRetries",
            input_type=Int(optional=True),
            prefix="-gcs-retries",
            doc="(--gcs-max-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20. ",
        ),
        ToolInput(
            tag="gcsProjectForRequesterPays",
            input_type=String(optional=True),
            prefix="--gcs-project-for-requester-pays",
            doc=" Project to bill when accessing requester pays  buckets. If unset, these buckets cannot be accessed.  Default value: . ",
        ),
        ToolInput(
            tag="intervalExclusionPadding",
            input_type=Int(optional=True),
            prefix="--interval-exclusion-padding",
            doc="(-ixp:Integer)  Amount of padding (in bp) to add to each interval you are excluding.  Default value: 0. ",
        ),
        ToolInput(
            tag="imr",
            input_type=String(optional=True),
            prefix="-imr:IntervalMergingRule",
            doc="(--interval-merging-rule)  Interval merging rule for abutting intervals  Default value: ALL. Possible values: {ALL, OVERLAPPING_ONLY} ",
        ),
        ToolInput(
            tag="ip",
            input_type=Int(optional=True),
            prefix="-ip",
            doc="(--interval-padding) Default value: 0.",
        ),
        ToolInput(
            tag="isr",
            input_type=String(optional=True),
            prefix="-isr:IntervalSetRule",
            doc="(--interval-set-rule)  Set merging approach to use for combining interval inputs  Default value: UNION. Possible values: {UNION, INTERSECTION} ",
        ),
        ToolInput(
            tag="le",
            input_type=Boolean(optional=True),
            prefix="--lenient",
            doc="(-LE) Lenient processing of VCF files Default value: false. Possible values: {true, false}",
        ),
        ToolInput(
            tag="quiet",
            input_type=Boolean(optional=True),
            prefix="--QUIET",
            doc="Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="readFilter",
            input_type=String(optional=True),
            prefix="--read-filter",
            doc="(-RF:String) Read filters to be applied before analysis This argument may be specified 0 or more times. Default value: null. Possible Values: {AlignmentAgreesWithHeaderReadFilter, AllowAllReadsReadFilter, AmbiguousBaseReadFilter, CigarContainsNoNOperator, FirstOfPairReadFilter, FragmentLengthReadFilter, GoodCigarReadFilter, HasReadGroupReadFilter, IntervalOverlapReadFilter, LibraryReadFilter, MappedReadFilter, MappingQualityAvailableReadFilter, MappingQualityNotZeroReadFilter, MappingQualityReadFilter, MatchingBasesAndQualsReadFilter, MateDifferentStrandReadFilter, MateOnSameContigOrNoMappedMateReadFilter, MateUnmappedAndUnmappedReadFilter, MetricsReadFilter, NonChimericOriginalAlignmentReadFilter, NonZeroFragmentLengthReadFilter, NonZeroReferenceLengthAlignmentReadFilter, NotDuplicateReadFilter, NotOpticalDuplicateReadFilter, NotSecondaryAlignmentReadFilter, NotSupplementaryAlignmentReadFilter, OverclippedReadFilter, PairedReadFilter, PassesVendorQualityCheckReadFilter, PlatformReadFilter, PlatformUnitReadFilter, PrimaryLineReadFilter, ProperlyPairedReadFilter, ReadGroupBlackListReadFilter, ReadGroupReadFilter, ReadLengthEqualsCigarLengthReadFilter, ReadLengthReadFilter, ReadNameReadFilter, ReadStrandFilter, SampleReadFilter, SecondOfPairReadFilter, SeqIsStoredReadFilter, SoftClippedReadFilter, ValidAlignmentEndReadFilter, ValidAlignmentStartReadFilter, WellformedReadFilter}",
        ),
        ToolInput(
            tag="readIndex",
            input_type=String(optional=True),
            prefix="-read-index",
            doc="(--read-index)  Indices to use for the read inputs. If specified, an index must be provided for every read input and in the same order as the read inputs. If this argument is not specified, the path to the index for each input will be inferred automatically.  This argument may be specified 0 or more times. Default value: null. ",
        ),
        ToolInput(
            tag="readValidationStringency",
            input_type=String(optional=True),
            prefix="--read-validation-stringency",
            doc="(-VS:ValidationStringency)  Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: SITool returned: 0 LENT. Possible values: {STRICT, LENIENT, SILENT} ",
        ),
        ToolInput(
            tag="reference",
            input_type=FastaWithDict(optional=True),
            prefix="--reference",
            doc="(-R:String) Reference sequence Default value: null.",
        ),
        ToolInput(
            tag="secondsBetweenProgressUpdates",
            input_type=Double(optional=True),
            prefix="-seconds-between-progress-updates",
            doc="(--seconds-between-progress-updates)  Output traversal statistics every time this many seconds elapse  Default value: 10.0. ",
        ),
        ToolInput(
            tag="sequenceDictionary",
            input_type=String(optional=True),
            prefix="-sequence-dictionary",
            doc="(--sequence-dictionary)  Use the given sequence dictionary as the master/canonical sequence dictionary.  Must be a .dict file.  Default value: null. ",
        ),
        ToolInput(
            tag="sitesOnlyVcfOutput",
            input_type=Boolean(optional=True),
            prefix="--sites-only-vcf-output:Boolean",
            doc=" If true, don't emit genotype fields when writing vcf file output.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="splitLibraryName",
            input_type=String(optional=True),
            prefix="--split-library-name",
            doc="(-LB)  Split file by library.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="rg",
            input_type=String(optional=True),
            prefix="--split-read-group",
            doc="(-RG:BooleanSplit) Default value: false. Possible values: {true, false}",
        ),
        ToolInput(
            tag="splitSample",
            input_type=String(optional=True),
            prefix="--split-sample",
            doc="(-SM:Boolean) Split file by sample. Default value: false. Possible values: {true, false}",
        ),
        ToolInput(
            tag="tmpDir",
            input_type=String(optional=True),
            prefix="--tmp-dir:GATKPathSpecifier",
            doc="Temp directory to use. Default value: null.",
        ),
        ToolInput(
            tag="jdkDeflater",
            input_type=Boolean(optional=True),
            prefix="-jdk-deflater",
            doc="(--use-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="jdkInflater",
            input_type=Boolean(optional=True),
            prefix="-jdk-inflater",
            doc="(--use-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="verbosity",
            input_type=String(optional=True),
            prefix="-verbosity:LogLevel",
            doc="(--verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG} ",
        ),
        ToolInput(
            tag="disableToolDefaultReadFilters",
            input_type=Boolean(optional=True),
            prefix="-disable-tool-default-read-filters",
            doc="(--disable-tool-default-read-filters)  Disable all tool default read filters (WARNING: many tools will not function correctly without their default read filters on)  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="ambigFilterBases",
            input_type=Int(optional=True),
            prefix="--ambig-filter-bases",
            doc="Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise, overrides threshold fraction.  Default value: null.  Cannot be used in conjuction with argument(s) maxAmbiguousBaseFraction",
        ),
        ToolInput(
            tag="ambigFilterFrac",
            input_type=Double(optional=True),
            prefix="--ambig-filter-frac",
            doc="Threshold fraction of ambiguous bases Default value: 0.05. Cannot be used in conjuction with argument(s) maxAmbiguousBases",
        ),
        ToolInput(
            tag="maxFragmentLength",
            input_type=Int(optional=True),
            prefix="--max-fragment-length",
            doc="Default value: 1000000.",
        ),
        ToolInput(
            tag="minFragmentLength",
            input_type=Int(optional=True),
            prefix="--min-fragment-length",
            doc="Default value: 0.",
        ),
        ToolInput(
            tag="keepIntervals",
            input_type=String(optional=True),
            prefix="--keep-intervals",
            doc='Valid only if "IntervalOverlapReadFilter" is specified: One or more genomic intervals to keep This argument must be specified at least once. Required. ',
        ),
        ToolInput(
            tag="library",
            input_type=String(optional=True),
            prefix="-library",
            doc='(--library) Valid only if "LibraryReadFilter" is specified: Name of the library to keep This argument must be specified at least once. Required.',
        ),
        ToolInput(
            tag="maximumMappingQuality",
            input_type=Int(optional=True),
            prefix="--maximum-mapping-quality",
            doc=" Maximum mapping quality to keep (inclusive)  Default value: null. ",
        ),
        ToolInput(
            tag="minimumMappingQuality",
            input_type=Int(optional=True),
            prefix="--minimum-mapping-quality",
            doc=" Minimum mapping quality to keep (inclusive)  Default value: 10. ",
        ),
        ToolInput(
            tag="dontRequireSoftClipsBothEnds",
            input_type=Boolean(optional=True),
            prefix="--dont-require-soft-clips-both-ends",
            doc=" Allow a read to be filtered out based on having only 1 soft-clipped block. By default, both ends must have a soft-clipped block, setting this flag requires only 1 soft-clipped block  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="filterTooShort",
            input_type=Int(optional=True),
            prefix="--filter-too-short",
            doc="Minimum number of aligned bases Default value: 30.",
        ),
        ToolInput(
            tag="platformFilterName",
            input_type=String(optional=True),
            prefix="--platform-filter-name:String",
            doc="This argument must be specified at least once. Required.",
        ),
        ToolInput(
            tag="blackListedLanes",
            input_type=String(optional=True),
            prefix="--black-listed-lanes:String",
            doc="Platform unit (PU) to filter out This argument must be specified at least once. Required.",
        ),
        ToolInput(
            tag="readGroupBlackList",
            input_type=String(optional=True),
            prefix="--read-group-black-list:StringThe",
            doc="This argument must be specified at least once. Required. ",
        ),
        ToolInput(
            tag="keepReadGroup",
            input_type=String(optional=True),
            prefix="--keep-read-group:String",
            doc="The name of the read group to keep Required.",
        ),
        ToolInput(
            tag="maxReadLength",
            input_type=Int(optional=True),
            prefix="--max-read-length",
            doc="Keep only reads with length at most equal to the specified value Required.",
        ),
        ToolInput(
            tag="minReadLength",
            input_type=Int(optional=True),
            prefix="--min-read-length",
            doc="Keep only reads with length at least equal to the specified value Default value: 1.",
        ),
        ToolInput(
            tag="readName",
            input_type=String(optional=True),
            prefix="--read-name:String",
            doc="Keep only reads with this read name Required.",
        ),
        ToolInput(
            tag="keepReverseStrandOnly",
            input_type=Boolean(optional=True),
            prefix="--keep-reverse-strand-only",
            doc=" Keep only reads on the reverse strand  Required. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="sample",
            input_type=String(optional=True),
            prefix="-sample:String",
            doc="(--sample) The name of the sample(s) to keep, filtering out all others This argument must be specified at least once. Required. ",
        ),
        ToolInput(
            tag="invertSoftClipRatioFilter",
            input_type=Boolean(optional=True),
            prefix="--invert-soft-clip-ratio-filter",
            doc=" Inverts the results from this filter, causing all variants that would pass to fail and visa-versa.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            tag="softClippedLeadingTrailingRatio",
            input_type=Double(optional=True),
            prefix="--soft-clipped-leading-trailing-ratio",
            doc=" Threshold ratio of soft clipped bases (leading / trailing the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumSoftClippedRatio",
        ),
        ToolInput(
            tag="softClippedRatioThreshold",
            input_type=Double(optional=True),
            prefix="--soft-clipped-ratio-threshold",
            doc=" Threshold ratio of soft clipped bases (anywhere in the cigar string) to total bases in read for read to be filtered.  Default value: null.  Cannot be used in conjuction with argument(s) minimumLeadingTrailingSoftClippedRatio",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.recalibrated.bam",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                    "javaOptions": ["-Xmx3G"],
                    "outputFilename": ".",
                },
                output=BamBai.basic_test(
                    "out",
                    2600900,
                    21300,
                    f"{remote_dir}/NA12878-BRCA1.split.flagstat",
                ),
            )
        ]


class Gatk4SplitReads_4_1_2(Gatk_4_1_2_0, Gatk4SplitReadsBase):
    pass


class Gatk4SplitReads_4_1_3(Gatk_4_1_3_0, Gatk4SplitReadsBase):
    pass


class Gatk4SplitReads_4_1_4(Gatk_4_1_4_0, Gatk4SplitReadsBase):
    pass


Gatk4SplitReadsLatest = Gatk4SplitReads_4_1_4



### Gatk4CollectInsertSizeMetrics ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 16,
            CaptureType.THIRTYX: 32,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class Gatk4CollectInsertSizeMetricsBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "CollectInsertSizeMetrics"

    def tool(self):
        return "Gatk4CollectInsertSizeMetrics"

    def friendly_name(self):
        return "GATK4: CollectInsertSizeMetrics"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Jiaan Yu"],
            dateCreated=date(2020, 2, 17),
            dateUpdated=date(2020, 2, 17),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "picard", "CollectInsertSizeMetrics"],
            documentationUrl="https://gatk.broadinstitute.org/hc/en-us/articles/360036715591-CollectInsertSizeMetrics-Picard-",
            documentation="Provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries",
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

    def inputs(self):
        return [
            *super(Gatk4CollectInsertSizeMetricsBase, self).inputs(),
            ToolInput(
                "bam",
                BamBai(optional=False),
                prefix="-I",
                doc="Input SAM or BAM file.  Required.",
                position=10,
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("bam", remove_file_extension=True),
                    extension=".txt",
                    suffix=".metrics",
                ),
                prefix="-O",
                doc="File to write the output to.  Required.",
            ),
            ToolInput(
                "outputHistogram",
                Filename(
                    prefix=InputSelector("bam", remove_file_extension=True),
                    extension=".pdf",
                    suffix=".histogram",
                ),
                prefix="-H",
                doc="File to write insert size Histogram chart to.  Required. ",
            ),
            *Gatk4CollectInsertSizeMetricsBase.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput("out", TextFile(), glob=InputSelector("outputFilename")),
            ToolOutput(
                "outHistogram",
                File(extension=".pdf"),
                glob=InputSelector("outputHistogram"),
            ),
        ]

    additional_args = [
        ToolInput(
            "argumentsFile",
            Array(File(), optional=True),
            prefix="--arguments_file",
            position=10,
            prefix_applies_to_all_elements=True,
            doc="read one or more arguments files and add them to the command line",
        ),
        ToolInput(
            "assumeSorted",
            Boolean(optional=True),
            prefix="--ASSUME_SORTED",
            position=11,
            doc="If true (default), then the sort order in the header file will be ignored.  Default value: true. Possible values: {true, false}",
        ),
        ToolInput(
            "deviations",
            Double(optional=True),
            prefix="--DEVIATIONS",
            position=11,
            doc="Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.  Default value: 10.0. ",
        ),
        ToolInput(
            "histogramWidth",
            Int(optional=True),
            prefix="--HISTOGRAM_WIDTH",
            position=11,
            doc="Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.  Default value: null. ",
        ),
        ToolInput(
            "includeDuplicates",
            Boolean(optional=True),
            prefix="--INCLUDE_DUPLICATES",
            position=11,
            doc="If true, also include reads marked as duplicates in the insert size histogram.  Default value: false. Possible values: {true, false} ",
        ),
        ToolInput(
            "metricAccumulationLevel",
            String(optional=True),
            prefix="--METRIC_ACCUMULATION_LEVEL",
            position=11,
            doc="The level(s) at  which to accumulate metrics.    This argument may be specified 0 or more times. Default value: [ALL_READS]. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} .",
        ),
        ToolInput(
            "minimumPCT",
            Float(optional=True),
            prefix="--MINIMUM_PCT",
            position=11,
            doc="When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this percentage of overall reads. (Range: 0 to 1).  Default value: 0.05.",
        ),
        ToolInput(
            "stopAfter",
            Int(optional=True),
            prefix="--STOP_AFTER",
            position=11,
            doc="Stop after  processing N reads, mainly for debugging.  Default value: 0. ",
        ),
        ToolInput(
            "version",
            Boolean(optional=True),
            prefix="--version",
            position=11,
            doc="display the version number for this tool Default value: false. Possible values: {true, false}",
        ),
        ToolInput(
            "showHidden",
            Boolean(optional=True),
            prefix="--showHidden",
            position=11,
            doc="display hidden  arguments  Default  value: false.  Possible values: {true, false} ",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        # The first 5 lines of the file include headers that change with every run (time, etc)
        with open(
            os.path.join(
                BioinformaticsTool.test_data_path(),
                "NA12878-BRCA1.markduped.metrics.txt",
            ),
            "r",
        ) as f:
            for i in range(5):
                next(f)
            expected_content = f.read()
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "javaOptions": ["-Xmx6G"],
                },
                output=TextFile.basic_test("out", 7260, expected_content, 905)
                + [
                    TTestExpectedOutput(
                        tag="outHistogram",
                        preprocessor=TTestPreprocessor.FileSize,
                        operator=operator.ge,
                        expected_value=15600,
                    ),
                ],
            )
        ]

class Gatk4CollectInsertSizeMetrics_4_0(Gatk_4_0_12, Gatk4CollectInsertSizeMetricsBase):
    pass


class Gatk4CollectInsertSizeMetrics_4_1_2(
    Gatk_4_1_2_0, Gatk4CollectInsertSizeMetricsBase
):
    pass


class Gatk4CollectInsertSizeMetrics_4_1_3(
    Gatk_4_1_3_0, Gatk4CollectInsertSizeMetricsBase
):
    pass


class Gatk4CollectInsertSizeMetrics_4_1_4(
    Gatk_4_1_4_0, Gatk4CollectInsertSizeMetricsBase
):
    pass


Gatk4CollectInsertSizeMetricsLatest = Gatk4CollectInsertSizeMetrics_4_1_3




### APPLY_BQSR ###



class Gatk4ApplyBqsrBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "ApplyBQSR"

    def friendly_name(self):
        return "GATK4: Apply base quality score recalibration"

    def tool(self):
        return "Gatk4ApplyBQSR"

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

    def inputs(self):
        return [
            *super(Gatk4ApplyBqsrBase, self).inputs(),
            ToolInput(
                "bam",
                BamBai(),
                prefix="-I",
                doc="The SAM/BAM/CRAM file containing reads.",
                secondaries_present_as={".bai": "^.bai"},
                position=10,
            ),
            ToolInput(
                "reference", FastaWithDict(), prefix="-R", doc="Reference sequence"
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("bam", remove_file_extension=True),
                    suffix=".recalibrated",
                    extension=".bam",
                ),
                prefix="-O",
                doc="Write output to this file",
            ),
            ToolInput(
                "recalFile",
                Tsv(optional=True),
                prefix="--bqsr-recal-file",
                doc="Input recalibration table for BQSR",
            ),
            ToolInput(
                "intervals",
                Bed(optional=True),
                prefix="--intervals",
                doc="-L (BASE) One or more genomic intervals over which to operate",
            ),
            ToolInput(
                "intervalStrings",
                Array(String, optional=True),
                prefix="--intervals",
                prefix_applies_to_all_elements=True,
                doc="-L (BASE) One or more genomic intervals over which to operate",
            ),
            *self.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai(),
                glob=InputSelector("outputFilename"),
                secondaries_present_as={".bai": "^.bai"},
            )
        ]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php",
            documentation="""
Apply base quality score recalibration: This tool performs the second pass in a two-stage 
process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
base qualities of the input reads based on the recalibration table produced by the 
BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
that affect the assignment of base quality scores by the sequencer. The first pass consists 
of calculating error empirically and finding patterns in how error varies with basecall 
features over all bases. The relevant observations are written to a recalibration table. 
The second pass consists of applying numerical corrections to each individual basecall 
based on the patterns identified in the first step (recorded in the recalibration table) 
and write out the recalibrated data to a new BAM or CRAM file.

- This tool replaces the use of PrintReads for the application of base quality score 
    recalibration as practiced in earlier versions of GATK (2.x and 3.x).
- You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
- Original qualities can be retained in the output file under the "OQ" tag if desired. 
    See the `--emit-original-quals` argument for details.
""".strip(),
        )

    additional_args = [
        # Put more detail in here from documentation
        ToolInput(
            "tmpDir",
            String(optional=True),
            prefix="--tmp-dir",
            position=11,
            default="/tmp/",
            doc="Temp directory to use.",
        )
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "recalFile": f"{remote_dir}/NA12878-BRCA1.markduped.table",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                },
                output=BamBai.basic_test(
                    "out",
                    2600000,
                    21000,
                    f"{remote_dir}/NA12878-BRCA1.recalibrated.flagstat",
                ),
            ),
            TTestCase(
                name="minimal",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "recalFile":f"{remote_dir}/NA12878-BRCA1.markduped.table",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                },
                output=self.minimal_test(),
            ),
        ]


class Gatk4ApplyBqsr_4_0(Gatk_4_0_12, Gatk4ApplyBqsrBase):
    pass


class Gatk4ApplyBqsr_4_1_2(Gatk_4_1_2_0, Gatk4ApplyBqsrBase):
    pass


class Gatk4ApplyBqsr_4_1_3(Gatk_4_1_3_0, Gatk4ApplyBqsrBase):
    pass


class Gatk4ApplyBqsr_4_1_4(Gatk_4_1_4_0, Gatk4ApplyBqsrBase):
    pass


Gatk4ApplyBqsrLatest = Gatk4ApplyBqsr_4_1_3



### BASE_RECALIBRATOR ###


class Gatk4BaseRecalibratorBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "BaseRecalibrator"

    def friendly_name(self):
        return "GATK4: Base Recalibrator"

    def tool(self):
        return "Gatk4BaseRecalibrator"

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 1

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 16

    def inputs(self):
        return [
            *super(Gatk4BaseRecalibratorBase, self).inputs(),
            *Gatk4BaseRecalibratorBase.additional_args,
            ToolInput(
                "bam",
                BamBai(),
                position=6,
                prefix="-I",
                doc="BAM/SAM/CRAM file containing reads",
                secondaries_present_as={".bai": "^.bai"},
            ),
            ToolInput(
                "knownSites",
                Array(VcfTabix()),
                prefix="--known-sites",
                position=28,
                prefix_applies_to_all_elements=True,
                doc="**One or more databases of known polymorphic sites used to exclude "
                "regions around known polymorphisms from analysis.** "
                "This algorithm treats every reference mismatch as an indication of error. However, real "
                "genetic variation is expected to mismatch the reference, so it is critical that a "
                "database of known polymorphic sites is given to the tool in order to skip over those sites. "
                "This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as "
                "this database. For users wishing to exclude an interval list of known variation simply "
                "use -XL my.interval.list to skip over processing those sites. Please note however "
                "that the statistics reported by the tool will not accurately reflected those sites "
                "skipped by the -XL argument.",
            ),
            ToolInput(
                "reference",
                FastaWithDict(),
                position=5,
                prefix="-R",
                doc="Reference sequence file",
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("bam", remove_file_extension=True),
                    extension=".table",
                ),
                position=8,
                prefix="-O",
                doc="**The output recalibration table filename to create.** "
                "After the header, data records occur one per line until the end of the file. The first "
                "several items on a line are the values of the individual covariates and will change "
                "depending on which covariates were specified at runtime. The last three items are the "
                "data- that is, number of observations for this combination of covariates, number of "
                "reference mismatches, and the raw empirical quality score calculated by phred-scaling "
                "the mismatch rate. Use '/dev/stdout' to print to standard out.",
            ),
            ToolInput(
                "intervals",
                Bed(optional=True),
                prefix="--intervals",
                doc="-L (BASE) One or more genomic intervals over which to operate",
            ),
            ToolInput(
                "intervalStrings",
                Array(String, optional=True),
                prefix="--intervals",
                prefix_applies_to_all_elements=True,
                doc="-L (BASE) One or more genomic intervals over which to operate",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Tsv(), glob=InputSelector("outputFilename"))]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "base recalibrator"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php",
            documentation="""
First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
Since there is a large amount of data one can then calculate an empirical probability of error given the 
particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
table (of the several covariate values, num observations, num mismatches, empirical quality score).  
""".strip(),
        )

    additional_args = [
        ToolInput(
            "tmpDir",
            String(optional=True),
            prefix="--tmp-dir",
            default="/tmp/",
            doc="Temp directory to use.",
        )
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "knownSites": [
                        f"{remote_dir}/Homo_sapiens_assembly38.known_indels.BRCA1.vcf.gz",
                        f"{remote_dir}/Homo_sapiens_assembly38.dbsnp138.BRCA1.vcf.gz",
                        f"{remote_dir}/Mills_and_1000G_gold_standard.indels.hg38.BRCA1.vcf.gz",
                        f"{remote_dir}/1000G_phase1.snps.high_confidence.hg38.BRCA1.vcf.gz",
                    ],
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                    "javaOptions": ["-Xmx12G"],
                },
                output=TextFile.basic_test(
                    "out", 1131758, "#:GATKReport.v1.1:5", 10376
                ),
            )
        ]
    
class Gatk4BaseRecalibrator_4_0(Gatk_4_0_12, Gatk4BaseRecalibratorBase):
    pass


class Gatk4BaseRecalibrator_4_1_2(Gatk_4_1_2_0, Gatk4BaseRecalibratorBase):
    pass


class Gatk4BaseRecalibrator_4_1_3(Gatk_4_1_3_0, Gatk4BaseRecalibratorBase):
    pass


class Gatk4BaseRecalibrator_4_1_4(Gatk_4_1_4_0, Gatk4BaseRecalibratorBase):
    pass


Gatk4BaseRecalibratorLatest = Gatk4BaseRecalibrator_4_1_4



### MERGESAMFILES ###


class Gatk4MergeSamFilesBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "MergeSamFiles"

    def tool(self):
        return "Gatk4MergeSamFiles"

    def friendly_name(self):
        return "GATK4: Merge SAM Files"

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

    def inputs(self):
        return [
            *super().inputs(),
            ToolInput(
                "bams",
                Array(BamBai()),
                prefix="-I",
                prefix_applies_to_all_elements=True,
                doc="The SAM/BAM file to sort.",
                position=10,
            ),
            ToolInput(
                "sampleName", String(optional=True), doc="Used for naming purposes only"
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("sampleName"),
                    suffix=".merged",
                    extension=".bam",
                ),
                position=10,
                prefix="-O",
                doc="SAM/BAM file to write merged result to",
            ),
            *self.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai(),
                glob=InputSelector("outputFilename"),
                secondaries_present_as={".bai": "^.bai"},
            )
        ]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "merge", "sam"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_MergeSamFiles.php",
            documentation="Merges multiple SAM/BAM files into one file",
        )

    additional_args = [
        ToolInput(
            "argumentsFile",
            Array(File(), optional=True),
            prefix="--arguments_file",
            position=10,
            doc="read one or more arguments files and add them to the command line",
        ),
        ToolInput(
            "assumeSorted",
            Boolean(optional=True),
            prefix="-AS",
            doc="If true, assume that the input files are in the same sort order as the requested "
            "output sort order, even if their headers say otherwise.",
        ),
        ToolInput(
            "comment",
            Array(String(), optional=True),
            prefix="-CO",
            doc="Comment(s) to include in the merged output file's header.",
        ),
        ToolInput(
            "mergeSequenceDictionaries",
            Boolean(optional=True),
            prefix="-MSD",
            doc="Merge the sequence dictionaries",
        ),
        ToolInput(
            "sortOrder",
            String(optional=True),
            prefix="-SO",
            position=10,
            doc="The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of "
            "the following values: [unsorted, queryname, coordinate, duplicate, unknown]",
        ),
        ToolInput(
            "useThreading",
            Boolean(optional=True),
            prefix="--USE_THREADING",
            doc="Option to create a background thread to encode, compress and write to disk the output file. "
            "The threaded version uses about 20% more CPU and decreases runtime by "
            "~20% when writing out a compressed BAM file.",
        ),
        ToolInput(
            "compressionLevel",
            Int(optional=True),
            prefix="--COMPRESSION_LEVEL",
            position=11,
            doc="Compression level for all compressed files created (e.g. BAM and GELI).",
        ),
        ToolInput(
            "createIndex",
            Boolean(optional=True),
            prefix="--CREATE_INDEX",
            position=11,
            doc="Whether to create a BAM index when writing a coordinate-sorted BAM file.",
        ),
        ToolInput(
            "createMd5File",
            Boolean(optional=True),
            prefix="--CREATE_MD5_FILE",
            position=11,
            doc="Whether to create an MD5 digest for any BAM or FASTQ files created.",
        ),
        ToolInput(
            "maxRecordsInRam",
            Int(optional=True),
            prefix="--MAX_RECORDS_IN_RAM",
            position=11,
            doc="When writing SAM files that need to be sorted, this will specify the number of "
            "records stored in RAM before spilling to disk. Increasing this number reduces "
            "the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.",
        ),
        ToolInput(
            "quiet",
            Boolean(optional=True),
            prefix="--QUIET",
            position=11,
            doc="Whether to suppress job-summary info on System.err.",
        ),
        ToolInput(
            "reference",
            FastaWithDict(optional=True),
            prefix="--reference",
            position=11,
            doc="Reference sequence file.",
        ),
        ToolInput(
            "tmpDir",
            String(optional=True),
            prefix="--TMP_DIR",
            position=11,
            default="/tmp/",
            doc="Undocumented option",
        ),
        ToolInput(
            "useJdkDeflater",
            Boolean(optional=True),
            prefix="--use_jdk_deflater",
            position=11,
            doc="Whether to use the JdkDeflater (as opposed to IntelDeflater)",
        ),
        ToolInput(
            "useJdkInflater",
            Boolean(optional=True),
            prefix="--use_jdk_inflater",
            position=11,
            doc="Whether to use the JdkInflater (as opposed to IntelInflater)",
        ),
        ToolInput(
            "validationStringency",
            String(optional=True),
            prefix="--VALIDATION_STRINGENCY",
            position=11,
            doc="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
            "can improve performance when processing a BAM file in which variable-length data "
            "(read, qualities, tags) do not otherwise need to be decoded."
            "The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), "
            "which can have one of the following values: [STRICT, LENIENT, SILENT]",
        ),
        ToolInput(
            "verbosity",
            String(optional=True),
            prefix="--verbosity",
            position=11,
            doc="The --verbosity argument is an enumerated type (LogLevel), which can have "
            "one of the following values: [ERROR, WARNING, INFO, DEBUG]",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bams": [
                        f"{remote_dir}/NA12878-BRCA1.sorted.bam",
                    ],
                    "createIndex": True,
                    "validationStringency": "SILENT",
                    "javaOptions": ["-Xmx6G"],
                    "maxRecordsInRam": 5000000,
                    "tmpDir": "./tmp",
                    "useThreading": True,
                },
                output=BamBai.basic_test(
                    "out",
                    2826968,
                    49688,
                    f"{remote_dir}/NA12878-BRCA1.bam.flagstat",
                    "963a51f7feed5b829319b947961b8a3e",
                    "231c10d0e43766170f5a7cd1b8a6d14e",
                ),
            )
        ]


class Gatk4MergeSamFiles_4_0(Gatk_4_0_12, Gatk4MergeSamFilesBase):
    pass


class Gatk4MergeSamFiles_4_1_2(Gatk_4_1_2_0, Gatk4MergeSamFilesBase):
    pass


class Gatk4MergeSamFiles_4_1_3(Gatk_4_1_3_0, Gatk4MergeSamFilesBase):
    pass


class Gatk4MergeSamFiles_4_1_4(Gatk_4_1_4_0, Gatk4MergeSamFilesBase):
    pass


Gatk4MergeSamFilesLatest = Gatk4MergeSamFiles_4_1_4


### MARKDUPLICATES ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 16,
            CaptureType.CHROMOSOME: 64,
            CaptureType.EXOME: 64,
            CaptureType.THIRTYX: 64,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class Gatk4MarkDuplicatesBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "MarkDuplicates"

    def tool(self):
        return "Gatk4MarkDuplicates"

    def friendly_name(self):
        return "GATK4: Mark Duplicates"

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

    def inputs(self):
        # Would be good to include this in the prefix:
        #   If(InputSelector("bam").length().equals(1), InputSelector("bam")[0].basename(), None)

        prefix = FirstOperator([InputSelector("outputPrefix"), "generated"])
        return [
            ToolInput(
                "bam",
                Array(Bam),
                prefix="-I",
                position=10,
                # secondaries_present_as={".bai": "^.bai"},
                doc="One or more input SAM or BAM files to analyze. Must be coordinate sorted.",
            ),
            ToolInput("outputPrefix", String(optional=True)),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=prefix,
                    suffix=".markduped",
                    extension=".bam",
                ),
                position=10,
                prefix="-O",
                doc="File to write duplication metrics to",
            ),
            ToolInput(
                "metricsFilename",
                Filename(prefix=prefix, suffix=".metrics", extension=".txt"),
                position=10,
                prefix="-M",
                doc="The output file to write marked records to.",
            ),
            *super().inputs(),
            *self.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai,
                glob=InputSelector("outputFilename"),
                secondaries_present_as={".bai": "^.bai"},
            ),
            ToolOutput("metrics", Tsv(), glob=InputSelector("metricsFilename")),
        ]

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "mark", "duplicates"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php",
            documentation="""MarkDuplicates (Picard): Identifies duplicate reads.

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are 
defined as originating from a single fragment of DNA. Duplicates can arise during sample 
preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for 
additional notes on PCR duplication artifacts. Duplicate reads can also result from a single 
amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the 
sequencing instrument. These duplication artifacts are referred to as optical duplicates.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads 
and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate
marking using molecular barcodes. After duplicate reads are collected, the tool differentiates 
the primary and duplicate reads using an algorithm that ranks reads by the sums of their 
base-quality scores (default method).

The tool's main output is a new SAM or BAM file, in which duplicates have been identified 
in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, 
which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, 
please see the following blog post for additional information.

Although the bitwise flag annotation indicates whether a read was marked as a duplicate, 
it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) 
tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. 
Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), 
only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the 
output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), 
as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 
This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the 
primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to 
skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are 
extremely large and estimating library complexity is not an aim. Note that without optical 
duplicate counts, library size estimation will be inaccurate.

MarkDuplicates also produces a metrics file indicating the numbers 
of duplicates for both single- and paired-end reads.

The program can take either coordinate-sorted or query-sorted inputs, however the behavior 
is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records 
and supplementary/secondary alignments are not marked as duplicates. However, when the input 
is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary 
reads are not excluded from the duplication test and can be marked as duplicate reads.

If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.""".strip(),
        )

    additional_args = [
        ToolInput(
            "argumentsFile",
            Array(File(), optional=True),
            prefix="--arguments_file",
            position=10,
            doc="read one or more arguments files and add them to the command line",
        ),
        ToolInput(
            "assumeSortOrder",
            String(optional=True),
            prefix="-ASO",
            doc="If not null, assume that the input file has this order even if the header says otherwise. "
            "Exclusion: This argument cannot be used at the same time as ASSUME_SORTED. "
            "The --ASSUME_SORT_ORDER argument is an enumerated type (SortOrder), which can have one of "
            "the following values: [unsorted, queryname, coordinate, duplicate, unknown]",
        ),
        ToolInput(
            "barcodeTag",
            String(optional=True),
            prefix="--BARCODE_TAG",
            doc="Barcode SAM tag (ex. BC for 10X Genomics)",
        ),
        ToolInput(
            "comment",
            Array(String(), optional=True),
            prefix="-CO",
            doc="Comment(s) to include in the output file's header.",
        ),
        # ToolInput(
        #     "compressionLevel",
        #     Int(optional=True),
        #     prefix="--COMPRESSION_LEVEL",
        #     position=11,
        #     doc="Compression level for all compressed files created (e.g. BAM and GELI).",
        # ),
        ToolInput(
            "createIndex",
            Boolean(optional=True),
            prefix="--CREATE_INDEX",
            default=True,
            position=11,
            doc="Whether to create a BAM index when writing a coordinate-sorted BAM file.",
        ),
        ToolInput(
            "createMd5File",
            Boolean(optional=True),
            prefix="--CREATE_MD5_FILE",
            position=11,
            doc="Whether to create an MD5 digest for any BAM or FASTQ files created.",
        ),
        ToolInput(
            "maxRecordsInRam",
            Int(optional=True),
            prefix="--MAX_RECORDS_IN_RAM",
            position=11,
            doc="When writing SAM files that need to be sorted, this will specify the number of "
            "records stored in RAM before spilling to disk. Increasing this number reduces "
            "the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.",
        ),
        ToolInput(
            "quiet",
            Boolean(optional=True),
            prefix="--QUIET",
            position=11,
            doc="Whether to suppress job-summary info on System.err.",
        ),
        ToolInput(
            "tmpDir",
            String(optional=True),
            prefix="--TMP_DIR",
            position=11,
            default="tmp/",
            doc="Undocumented option",
        ),
        ToolInput(
            "useJdkDeflater",
            Boolean(optional=True),
            prefix="--use_jdk_deflater",
            position=11,
            doc="Whether to use the JdkDeflater (as opposed to IntelDeflater)",
        ),
        ToolInput(
            "useJdkInflater",
            Boolean(optional=True),
            prefix="--use_jdk_inflater",
            position=11,
            doc="Whether to use the JdkInflater (as opposed to IntelInflater)",
        ),
        ToolInput(
            "validationStringency",
            String(optional=True),
            prefix="--VALIDATION_STRINGENCY",
            position=11,
            doc="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
            "can improve performance when processing a BAM file in which variable-length data "
            "(read, qualities, tags) do not otherwise need to be decoded."
            "The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), "
            "which can have one of the following values: [STRICT, LENIENT, SILENT]",
        ),
        ToolInput(
            "verbosity",
            String(optional=True),
            prefix="--verbosity",
            position=11,
            doc="The --verbosity argument is an enumerated type (LogLevel), which can have "
            "one of the following values: [ERROR, WARNING, INFO, DEBUG]",
        ),
        ToolInput(
            "opticalDuplicatePixelDistance",
            Int(optional=True),
            prefix="--OPTICAL_DUPLICATE_PIXEL_DISTANCE",
            doc="The maximum offset between two duplicate clusters in order to consider them optical duplicates. "
            "The default is appropriate for unpatterned versions of the Illumina platform. For the patterned "
            "flowcell models, 2500 is more appropriate. For other platforms and models, users should experiment "
            "to find what works best.",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": [
                        f"{remote_dir}/NA12878-BRCA1.merged.bam",
                    ],
                    "javaOptions": ["-Xmx6G"],
                    "maxRecordsInRam": 5000000,
                    "createIndex": True,
                    "tmpDir": "./tmp",
                },
                output=BamBai.basic_test(
                    "out",
                    2829000,
                    3780,
                    f"{remote_dir}/NA12878-BRCA1.markduped.bam.flagstat",
                )
                + TextFile.basic_test(
                    "metrics",
                    3700,
                    "NA12878-BRCA1\t193\t9468\t164\t193\t46\t7\t1\t0.003137\t7465518",
                    112,
                ),
            )
        ]


class Gatk4MarkDuplicates_4_0(Gatk_4_0_12, Gatk4MarkDuplicatesBase):
    pass


class Gatk4MarkDuplicates_4_1_2(Gatk_4_1_2_0, Gatk4MarkDuplicatesBase):
    pass


class Gatk4MarkDuplicates_4_1_3(Gatk_4_1_3_0, Gatk4MarkDuplicatesBase):
    pass


class Gatk4MarkDuplicates_4_1_4(Gatk_4_1_4_0, Gatk4MarkDuplicatesBase):
    pass


Gatk4MarkDuplicatesLatest = Gatk4MarkDuplicates_4_1_4


### SORTSAM ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 16,
            CaptureType.THIRTYX: 32,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class Gatk4SortSamBase(Gatk4ToolBase, ABC):
    @classmethod
    def gatk_command(cls):
        return "SortSam"

    def tool(self):
        return "Gatk4SortSam"

    def friendly_name(self):
        return "GATK4: SortSAM"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="Broad Institute",
            doi=None,
            citation="See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information",
            keywords=["gatk", "gatk4", "broad", "sort", "sam"],
            documentationUrl="https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.3/org_broadinstitute_hellbender_tools_picard_sam_SortSam.php",
            documentation="Sorts a SAM/BAM/CRAM file.",
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

    def inputs(self):
        return [
            *super(Gatk4SortSamBase, self).inputs(),
            ToolInput(
                "bam",
                Bam(),
                prefix="-I",
                doc="The SAM/BAM/CRAM file to sort.",
                position=10,
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("bam", remove_file_extension=True),
                    suffix=".sorted",
                    extension=".bam",
                ),
                position=10,
                prefix="-O",
                doc="The sorted SAM/BAM/CRAM output file.",
            ),
            ToolInput(
                "sortOrder",
                String(),
                prefix="-SO",
                position=10,
                doc="The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of "
                "the following values: [unsorted, queryname, coordinate, duplicate, unknown]",
            ),
            *Gatk4SortSamBase.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai(),
                glob=InputSelector("outputFilename"),
                secondaries_present_as={".bai": "^.bai"},
            )
        ]

    additional_args = [
        ToolInput(
            "argumentsFile",
            Array(File(), optional=True),
            prefix="--arguments_file",
            position=10,
            prefix_applies_to_all_elements=True,
            doc="read one or more arguments files and add them to the command line",
        ),
        ToolInput(
            "compressionLevel",
            Int(optional=True),
            prefix="--COMPRESSION_LEVEL",
            position=11,
            doc="Compression level for all compressed files created (e.g. BAM and GELI).",
        ),
        ToolInput(
            "createIndex",
            Boolean(optional=True),
            prefix="--CREATE_INDEX",
            default=True,
            position=11,
            doc="Whether to create a BAM index when writing a coordinate-sorted BAM file.",
        ),
        ToolInput(
            "createMd5File",
            Boolean(optional=True),
            prefix="--CREATE_MD5_FILE",
            position=11,
            doc="Whether to create an MD5 digest for any BAM or FASTQ files created.",
        ),
        ToolInput(
            "maxRecordsInRam",
            Int(optional=True),
            prefix="--MAX_RECORDS_IN_RAM",
            position=11,
            doc="When writing SAM files that need to be sorted, this will specify the number of "
            "records stored in RAM before spilling to disk. Increasing this number reduces "
            "the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.",
        ),
        ToolInput(
            "quiet",
            Boolean(optional=True),
            prefix="--QUIET",
            position=11,
            doc="Whether to suppress job-summary info on System.err.",
        ),
        ToolInput(
            "reference",
            FastaWithDict(optional=True),
            prefix="--reference",
            position=11,
            doc="Reference sequence file.",
        ),
        ToolInput(
            "tmpDir",
            String(optional=True),
            prefix="--TMP_DIR",
            position=11,
            doc="Undocumented option",
            default="/tmp/",
        ),
        ToolInput(
            "useJdkDeflater",
            Boolean(optional=True),
            prefix="--use_jdk_deflater",
            position=11,
            doc="Whether to use the JdkDeflater (as opposed to IntelDeflater)",
        ),
        ToolInput(
            "useJdkInflater",
            Boolean(optional=True),
            prefix="--use_jdk_inflater",
            position=11,
            doc="Whether to use the JdkInflater (as opposed to IntelInflater)",
        ),
        ToolInput(
            "validationStringency",
            String(optional=True),
            prefix="--VALIDATION_STRINGENCY",
            position=11,
            doc="Validation stringency for all SAM files read by this program. Setting stringency to SILENT "
            "can improve performance when processing a BAM file in which variable-length data "
            "(read, qualities, tags) do not otherwise need to be decoded."
            "The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), "
            "which can have one of the following values: [STRICT, LENIENT, SILENT]",
        ),
        ToolInput(
            "verbosity",
            String(optional=True),
            prefix="--verbosity",
            position=11,
            doc="The --verbosity argument is an enumerated type (LogLevel), which can have "
            "one of the following values: [ERROR, WARNING, INFO, DEBUG]",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.bam",
                    "sortOrder": "coordinate",
                    "createIndex": True,
                    "maxRecordsInRam": 5000000,
                    "tmpDir": "./tmp",
                    "validationStringency": "SILENT",
                    "javaOptions": ["-Xmx6G"],
                },
                output=BamBai.basic_test(
                    "out",
                    2826980,
                    49688,
                    f"{remote_dir}/NA12878-BRCA1.bam.flagstat",
                    "15eb0f8168b42e8ce3ab8b9bc9199e3c",
                    "a9042025f29f7a08e5f56ce8d11469a1",
                ),
            ),
            TTestCase(
                name="minimal",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.bam",
                    "sortOrder": "coordinate",
                    "createIndex": True,
                    "maxRecordsInRam": 5000000,
                    "tmpDir": "./tmp",
                    "validationStringency": "SILENT",
                    "javaOptions": ["-Xmx6G"],
                },
                output=self.minimal_test(),
            ),
        ]


class Gatk4SortSam_4_0(Gatk_4_0_12, Gatk4SortSamBase):
    pass


class Gatk4SortSam_4_1_2(Gatk_4_1_2_0, Gatk4SortSamBase):
    pass


class Gatk4SortSam_4_1_3(Gatk_4_1_3_0, Gatk4SortSamBase):
    pass


class Gatk4SortSam_4_1_4(Gatk_4_1_4_0, Gatk4SortSamBase):
    pass


Gatk4SortSamLatest = Gatk4SortSam_4_1_4

if __name__ == "__main__":
    print(Gatk4SortSam_4_0().help())
