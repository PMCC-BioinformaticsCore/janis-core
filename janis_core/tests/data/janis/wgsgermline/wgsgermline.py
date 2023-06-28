from datetime import date
from typing import Optional, List

from janis_core import (
    Array,
    File,
    Float,
    InputQualityType,
    Int,
    String,
    StringFormatter,
    WorkflowMetadata,
    InputDocumentation,
)
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)
from janis_core.operators.standard import FirstOperator

from janis_bioinformatics.data_types import (
    FastaWithDict,
    FastqGzPair,
    Bam,
    BamBai,
    Bed,
    BedTabix,
    Vcf,
    VcfTabix,
    CompressedVcf,
)
from janis_unix.data_types import TextFile, ZipFile
from janis_bioinformatics.tools.bioinformaticstoolbase import BioinformaticsWorkflow

from janis_unix.tools import UncompressArchive
from janis_bioinformatics.tools.babrahambioinformatics import FastQC_0_11_8
from janis_bioinformatics.tools.bcftools import BcfToolsSort_1_9, BcfToolsConcat_1_9
from janis_bioinformatics.tools.common import BwaAligner, MergeAndMarkBams_4_1_3
from janis_bioinformatics.tools.common import GATKBaseRecalBQSRWorkflow_4_1_3
from janis_bioinformatics.tools.htslib import BGZipLatest
from janis_bioinformatics.tools.papenfuss import Gridss_2_6_2
from janis_bioinformatics.tools.pmac import (
    AddBamStatsGermline_0_1_0,
    CombineVariants_0_0_8,
    GenerateGenomeFileForBedtoolsCoverage,
    GenerateIntervalsByChromosome,
    GenerateMantaConfig,
    GenerateVardictHeaderLines,
    ParseFastqcAdapters,
    PerformanceSummaryGenome_0_1_0,
)
from janis_bioinformatics.tools.variantcallers import (
    GatkGermlineVariantCaller_4_1_3,
    IlluminaGermlineVariantCaller,
    VardictGermlineVariantCaller,
)

from janis_pipelines.reference import WGS_INPUTS

INPUT_DOCS = {
    **WGS_INPUTS,
    "fastqs": {
        "doc": "An array of FastqGz pairs. These are aligned separately and merged "
        "to create higher depth coverages from multiple sets of reads",
        "quality": InputQualityType.user,
        "example": [
            ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"],
            ["sample1_R1-TOPUP.fastq.gz", "sample1_R2-TOPUP.fastq.gz"],
        ],
    },
    "sample_name": {
        "doc": "Sample name from which to generate the readGroupHeaderLine for BwaMem",
        "quality": InputQualityType.user,
        "example": "NA12878",
    },
    "bam": {
        "doc": "Input indexed bam (+ .bam.bai) to process. You only specify the primary sample.bam, and the index (eg: NA12878.bam.bai) will be picked up automatically.",
        "quality": InputQualityType.user,
        "example": "NA12878.bam",
    },
}


class WGSGermlineMultiCallers(BioinformaticsWorkflow):
    def id(self):
        return "WGSGermlineMultiCallers"

    def friendly_name(self):
        return "WGS Germline (Multi callers)"

    def version(self):
        return "1.4.0"

    ### PIPELINE CONSTRUCTOR
    def constructor(self):
        self.add_inputs()

        self.add_fastqc()
        self.add_trim_and_align_fastq()
        self.add_merge_and_markdups_bam()
        self.add_bam_qc(bam_source=self.merge_and_markdups.out)
        self.add_gridss(bam_source=self.merge_and_markdups.out)
        self.add_gatk_variantcaller(bam_source=self.merge_and_markdups.out)
        self.add_strelka_variantcaller(bam_source=self.merge_and_markdups.out)
        self.add_vardict_variantcaller(bam_source=self.merge_and_markdups.out)
        self.add_combine_variants(bam_source=self.merge_and_markdups.out)

    ### INPUTS
    def add_inputs(self):
        self.input("sample_name", String, doc=INPUT_DOCS["sample_name"])
        self.input("fastqs", Array(FastqGzPair), doc=INPUT_DOCS["fastqs"])

        self.add_inputs_for_reference()
        self.add_inputs_for_adapter_trimming()
        self.add_inputs_for_intervals()
        self.add_inputs_for_configuration()

    def add_inputs_for_reference(self):
        self.input("reference", FastaWithDict, doc=INPUT_DOCS["reference"])
        self.input("snps_dbsnp", VcfTabix, doc=INPUT_DOCS["snps_dbsnp"])
        self.input("snps_1000gp", VcfTabix, doc=INPUT_DOCS["snps_1000gp"])
        self.input("known_indels", VcfTabix, doc=INPUT_DOCS["known_indels"])
        self.input("mills_indels", VcfTabix, doc=INPUT_DOCS["mills_indels"])

    def add_inputs_for_adapter_trimming(self):
        self.input("adapter_file", File)
        self.input("contaminant_file", File)

    def add_inputs_for_intervals(self):
        self.input(
            "gatk_intervals",
            Array(Bed, optional=True),
            doc=INPUT_DOCS["gatk_intervals"],
        )
        self.input("vardict_intervals", Array(Bed), doc=INPUT_DOCS["vardict_intervals"])
        self.input("strelka_intervals", BedTabix, doc=INPUT_DOCS["strelka_intervals"])
        self.input("gridss_blacklist", Bed, doc=INPUT_DOCS["gridss_blacklist"])

    def add_inputs_for_configuration(self):
        # vardict related options
        self.input(
            "allele_freq_threshold",
            Float,
            0.05,
            doc=InputDocumentation(
                "The threshold for VarDict's allele frequency, default: 0.05 or 5%",
                quality=InputQualityType.configuration,
            ),
        ),
        self.input("minMappingQual", Int(optional=True))
        self.input("filter", String(optional=True))

    ### PIPELINE STEPS
    def add_fastqc(self):
        self.step(
            "fastqc",
            FastQC_0_11_8(reads=self.fastqs),
            scatter="reads",
        )
        self.output(
            "out_fastqc_R1_reports",
            source=self.fastqc.out_R1,
            output_folder="reports",
            doc="A zip file of the FastQC quality report.",
        )
        self.output(
            "out_fastqc_R2_reports",
            source=self.fastqc.out_R2,
            output_folder="reports",
            doc="A zip file of the FastQC quality report.",
        )

    def add_trim_and_align_fastq(self):
        self.step(
            "getfastqc_adapters",
            ParseFastqcAdapters(
                read1_fastqc_datafile=self.fastqc.out_R1_datafile,
                read2_fastqc_datafile=self.fastqc.out_R2_datafile,
                adapters_lookup=self.adapter_file,
                contamination_lookup=self.contaminant_file,
            ),
            scatter=["read1_fastqc_datafile", "read2_fastqc_datafile"],
        )
        self.step(
            "align_and_sort",
            BwaAligner(
                fastq=self.fastqs,
                reference=self.reference,
                sample_name=self.sample_name,
                sortsam_tmpDir="./tmp",
                three_prime_adapter_read1=self.getfastqc_adapters.out_R1_sequences,
                three_prime_adapter_read2=self.getfastqc_adapters.out_R2_sequences,
            ),
            scatter=["fastq", "three_prime_adapter_read1", "three_prime_adapter_read2"],
        )

    def add_merge_and_markdups_bam(self):
        self.step(
            "merge_and_markdups",
            MergeAndMarkBams_4_1_3(
                bams=self.align_and_sort.out, sampleName=self.sample_name
            ),
        )
        self.output(
            "out_bam",
            source=self.merge_and_markdups.out,
            output_folder="bams",
            doc="Aligned and indexed bam.",
            output_name=self.sample_name,
        )

    def add_bam_qc(self, bam_source):
        # Temporarily remove GATK4 DepthOfCoverage for performance reasons, see:
        #   https://gatk.broadinstitute.org/hc/en-us/community/posts/360071895391-Speeding-up-GATK4-DepthOfCoverage
        # Gatk4DepthOfCoverage
        # self.step(
        #     "coverage",
        #     Gatk4DepthOfCoverage_4_1_6(
        #         bam=bam_source,
        #         reference=self.reference,
        #         outputPrefix=self.sample_name,
        #         intervals=intervals,
        #         # current version gatk 4.1.6.0 only support --count-type as COUNT_READS
        #         # countType="COUNT_FRAGMENTS_REQUIRE_SAME_BASE",
        #         omitDepthOutputAtEachBase=True,
        #         summaryCoverageThreshold=[1, 50, 100, 300, 500],
        #     ),
        # )
        self.step(
            "calculate_performancesummary_genomefile",
            GenerateGenomeFileForBedtoolsCoverage(reference=self.reference),
        )
        self.step(
            "performance_summary",
            PerformanceSummaryGenome_0_1_0(
                bam=bam_source,
                genome_file=self.calculate_performancesummary_genomefile.out,
                sample_name=self.sample_name,
            ),
        )
        # COVERGAE
        # self.output(
        #     "sample_coverage",
        #     source=self.coverage.out_sampleSummary,
        #     output_folder=["performance_summary", self.sample_name],
        #     doc="A text file of depth of coverage summary of bam",
        # )
        # BAM PERFORMANCE
        self.output(
            "out_performance_summary",
            source=self.performance_summary.performanceSummaryOut,
            output_folder=["performance_summary", self.sample_name],
            doc="A text file of performance summary of bam",
        )

    def add_gridss(self, bam_source):
        # GRIDSS
        self.step(
            "vc_gridss",
            Gridss_2_6_2(
                bams=[bam_source],
                reference=self.reference,
                blacklist=self.gridss_blacklist,
            ),
        )
        self.output(
            "out_gridss_assembly",
            source=self.vc_gridss.assembly,
            output_folder="gridss",
            doc="Assembly returned by GRIDSS",
        )
        self.output(
            "out_variants_gridss",
            source=self.vc_gridss.out,
            output_folder="gridss",
            doc="Variants from the GRIDSS variant caller",
        )

    def add_gatk_variantcaller(self, bam_source):
        # CREATE GATK INTERVAL
        intervals = FirstOperator(
            [
                self.gatk_intervals,
                self.step(
                    "generate_gatk_intervals",
                    GenerateIntervalsByChromosome(reference=self.reference),
                    when=self.gatk_intervals.is_null(),
                ).out_regions,
            ]
        )
        self.step(
            "bqsr",
            GATKBaseRecalBQSRWorkflow_4_1_3(
                bam=bam_source,
                reference=self.reference,
                snps_dbsnp=self.snps_dbsnp,
                snps_1000gp=self.snps_1000gp,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                intervals=intervals,
            ),
            scatter="intervals",
            doc="Perform base quality score recalibration",
        )
        self.step(
            "vc_gatk",
            GatkGermlineVariantCaller_4_1_3(
                bam=self.bqsr.out,
                intervals=intervals,
                reference=self.reference,
                snps_dbsnp=self.snps_dbsnp,
            ),
            scatter=["intervals", "bam"],
        )
        # POST VARIANT CALLING
        self.step(
            "vc_gatk_merge",
            BcfToolsConcat_1_9(vcf=self.vc_gatk.out.as_type(Array(Vcf))),
        )
        self.step(
            "vc_gatk_sort_combined",
            BcfToolsSort_1_9(vcf=self.vc_gatk_merge.out.as_type(CompressedVcf)),
        )
        self.step(
            "vc_gatk_uncompress",
            UncompressArchive(
                file=self.vc_gatk_sort_combined.out,
                force=True
            ),
        )
        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_sort_combined.out,
            output_folder=[
                "variants",
            ],
            output_name=StringFormatter(
                "{sample_name}_gatk",
                sample_name=self.sample_name,
            ),
            doc="Merged variants from the GATK caller",
        )
        self.output(
            "out_variants_gatk_split",
            source=self.vc_gatk.out,
            output_folder=[
                "variants",
                "GatkByInterval",
            ],
            doc="Unmerged variants from the GATK caller (by interval)",
        )

    def add_strelka_variantcaller(self, bam_source):
        # STRELKA
        # MANTA CONFIG FOR WGS CALLING
        self.step("generate_manta_config", GenerateMantaConfig())
        self.step(
            "vc_strelka",
            IlluminaGermlineVariantCaller(
                bam=bam_source,
                reference=self.reference,
                intervals=self.strelka_intervals,
                manta_config=self.generate_manta_config.out,
            ),
        )
        # POST VARIANT CALLING
        self.step("vc_strelka_compress", BGZipLatest(file=self.vc_strelka.out))
        self.output(
            "out_variants_strelka",
            source=self.vc_strelka_compress.out.as_type(CompressedVcf),
            output_folder=[
                "variants",
            ],
            output_name=StringFormatter(
                "{sample_name}_strelka",
                sample_name=self.sample_name,
            ),
            doc="Variants from the Strelka variant caller",
        )

    def add_vardict_variantcaller(self, bam_source):
        # VARDICT
        self.step(
            "generate_vardict_headerlines",
            GenerateVardictHeaderLines(reference=self.reference),
        )
        self.step(
            "vc_vardict",
            VardictGermlineVariantCaller(
                bam=bam_source,
                reference=self.reference,
                intervals=self.vardict_intervals,
                sample_name=self.sample_name,
                allele_freq_threshold=self.allele_freq_threshold,
                header_lines=self.generate_vardict_headerlines.out,
                minMappingQual=self.minMappingQual,
                filter=self.filter,
            ),
            scatter="intervals",
        )
        # POST VARIANT CALLING
        self.step(
            "vc_vardict_merge",
            BcfToolsConcat_1_9(vcf=self.vc_vardict.out.as_type(Array(Vcf))),
        )
        self.step(
            "vc_vardict_sort_combined",
            BcfToolsSort_1_9(vcf=self.vc_vardict_merge.out.as_type(CompressedVcf)),
        )
        self.step(
            "vc_vardict_uncompress_for_combine",
            UncompressArchive(
                file=self.vc_vardict_sort_combined.out,
                force=True
            ),
        )
        self.output(
            "out_variants_vardict",
            source=self.vc_vardict_sort_combined.out,
            output_folder=[
                "variants",
            ],
            output_name=StringFormatter(
                "{sample_name}_vardict",
                sample_name=self.sample_name,
            ),
            doc="Merged variants from the VarDict caller",
        )
        self.output(
            "out_variants_vardict_split",
            source=self.vc_vardict.out,
            output_folder=[
                "variants",
                "VardictByInterval",
            ],
            doc="Unmerged variants from the VarDict caller (by interval)",
        )

    def add_combine_variants(self, bam_source):
        # COMBINE VARIANTS FROM GATK, STRELKA AND VARDICT
        # Note, this is reliant on the specific step names from previous steps
        self.step(
            "combine_variants",
            CombineVariants_0_0_8(
                vcfs=[
                    self.vc_gatk_uncompress.out.as_type(Vcf),
                    self.vc_strelka.out,
                    self.vc_vardict_uncompress_for_combine.out.as_type(Vcf),
                ],
                type="germline",
                columns=["AC", "AN", "AF", "AD", "DP", "GT"],
            ),
        )
        self.step("combined_compress", BGZipLatest(file=self.combine_variants.out))
        self.step(
            "combined_sort",
            BcfToolsSort_1_9(vcf=self.combined_compress.out.as_type(CompressedVcf)),
        )
        self.step(
            "combined_uncompress", 
            UncompressArchive(
                file=self.combined_sort.out,
                force=True
            )
        )

        self.step(
            "combined_addbamstats",
            AddBamStatsGermline_0_1_0(
                bam=bam_source,
                vcf=self.combined_uncompress.out.as_type(Vcf),
                reference=self.reference,
            ),
        )
        self.output(
            "out_variants_combined_bamstats",
            source=self.combined_addbamstats.out,
            output_folder="variants",
            output_name=StringFormatter(
                "{sample_name}_combined",
                sample_name=self.sample_name,
            ),
            doc="Combined variants from all 3 callers",
        )

    ### UNIT TEST
    def tests(self) -> Optional[List[TTestCase]]:
        parent_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics"
        brca1_test_data = f"{parent_dir}/brca1_test/test_data"

        return [
            TTestCase(
                name="brca1",
                input={
                    "sample_name": "NA12878-BRCA1",
                    "fastqs": [
                        [
                            f"{brca1_test_data}/NA12878-BRCA1_R1.fastq.gz",
                            f"{brca1_test_data}/NA12878-BRCA1_R2.fastq.gz",
                        ]
                    ],
                    "reference": f"{brca1_test_data}/Homo_sapiens_assembly38.chr17.fasta",
                    "gridss_blacklist": f"{brca1_test_data}/consensusBlacklist.hg38.chr17.bed",
                    "gatk_intervals": [f"{brca1_test_data}/BRCA1.hg38.bed"],
                    "strelka_intervals": f"{brca1_test_data}/BRCA1.hg38.bed.gz",
                    "vardict_intervals": [
                        f"{brca1_test_data}/BRCA1.hg38.split-intervals.bed"
                    ],
                    "known_indels": f"{brca1_test_data}/Homo_sapiens_assembly38.known_indels.BRCA1.vcf.gz",
                    "mills_indels": f"{brca1_test_data}/Mills_and_1000G_gold_standard.indels.hg38.BRCA1.vcf.gz",
                    "snps_1000gp": f"{brca1_test_data}/1000G_phase1.snps.high_confidence.hg38.BRCA1.vcf.gz",
                    "snps_dbsnp": f"{brca1_test_data}/Homo_sapiens_assembly38.dbsnp138.BRCA1.vcf.gz",
                    "contaminant_file": f"{brca1_test_data}/contaminant_list.txt",
                    "adapter_file": f"{brca1_test_data}/adapter_list.txt",
                    "vc_strelka_manta_runtime_cpu": 1,
                    "vc_strelka_strelka_runtime_cpu": 1,
                    "vc_strelka_manta_runtime_memory": 5,
                    "vc_strelka_strelka_runtime_memory": 3,
                },
                output=Array.array_wrapper(
                    [ZipFile.basic_test("out_fastqc_R1_reports", 408000)]
                )
                + Array.array_wrapper(
                    [ZipFile.basic_test("out_fastqc_R2_reports", 408000)]
                )
                + BamBai.basic_test("out_bam", 2822000, 49600)
                + TextFile.basic_test(
                    "out_performance_summary",
                    948,
                    md5="575354942cfb8d0367725f9020181443",
                )
                + Bam.basic_test("out_gridss_assembly", 16000)
                + Vcf.basic_test("out_variants_gridss", 32000, 130)
                + CompressedVcf.basic_test("out_variants_gatk", 11000, 223)
                + Array.array_wrapper(
                    [Vcf.basic_test("out_variants_gatk_split", 51000, 221)]
                )
                + CompressedVcf.basic_test("out_variants_strelka", 7500, 226)
                + CompressedVcf.basic_test("out_variants_vardict", 19000, 260)
                + Array.array_wrapper(
                    [Vcf.basic_test("out_variants_vardict_split", 84000, 258)]
                )
                + Vcf.basic_test("out_variants_combined_bamstats", 99000, 307),
            )
        ]

    def bind_metadata(self):
        meta: WorkflowMetadata = self.metadata

        meta.keywords = [
            "wgs",
            "cancer",
            "germline",
            "variants",
            "gatk",
            "strelka",
            "vardict",
            "gridss",
        ]
        meta.contributors = ["Richard Lupat", "Michael Franklin", "Jiaan Yu"]
        meta.dateCreated = date(2018, 12, 24)
        meta.dateUpdated = date(2022, 2, 25)
        meta.short_documentation = (
            "A variant-calling WGS pipeline using GATK, VarDict and Strelka2."
        )
        meta.documentation = """\
This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs and call variants using:
This workflow is a reference pipeline using the Janis Python framework (pipelines assistant).
- Takes raw sequence data in the FASTQ format;
- Align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate variant callers (GRIDSS / GATK / Strelka / VarDict);
- Merges the variants from GATK / Strelka / VarDict.
- Outputs the final variants in the VCF format.
**Resources**
This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:
- https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/
This pipeline expects the assembly references to be as they appear in that storage \
    (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.
"""


if __name__ == "__main__":
    # import os.path

    # w = WGSGermlineMultiCallers()
    # args = {
    #     "to_console": False,
    #     "to_disk": False,
    #     "validate": True,
    #     "export_path": os.path.join(
    #         os.path.dirname(os.path.realpath(__file__)), "{language}"
    #     ),
    #     "with_resource_overrides": True,
    # }
    # w.translate("cwl", **args)
    # w.translate("wdl", **args)
    WGSGermlineMultiCallers().translate("wdl")
