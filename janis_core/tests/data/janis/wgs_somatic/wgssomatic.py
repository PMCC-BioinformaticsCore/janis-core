from typing import Optional, List
from datetime import date
from janis_core import (
    Array,
    Boolean,
    File,
    InputQualityType,
    Int,
    String,
    StringFormatter,
    WorkflowMetadata,
)
from janis_core.operators.standard import FirstOperator
from janis_core.tool.test_classes import TTestCase

from janis_unix.data_types import TextFile, ZipFile
from janis_bioinformatics.data_types import (
    FastqGzPair,
    Bam,
    BamBai,
    Vcf,
    CompressedVcf,
    VcfTabix,
)

from janis_unix.tools import UncompressArchive
from janis_bioinformatics.tools.bcftools import BcfToolsSort_1_9, BcfToolsConcat_1_9
from janis_bioinformatics.tools.common import (
    GATKBaseRecalBQSRWorkflow_4_1_3,
    FacetsWorkflow,
)
from janis_bioinformatics.tools.htslib import BGZipLatest
from janis_bioinformatics.tools.papenfuss import Gridss_2_6_2
from janis_bioinformatics.tools.pmac import (
    AddBamStatsSomatic_0_1_0,
    CircosPlot_0_1_2,
    CombineVariants_0_0_8,
    GenerateGenomeFileForBedtoolsCoverage,
    GenerateIntervalsByChromosome,
    GenerateMantaConfig,
    GenerateVardictHeaderLines,
    PerformanceSummaryGenome_0_1_0,
)
from janis_bioinformatics.tools.variantcallers.illuminasomatic_strelka import (
    IlluminaSomaticVariantCaller,
)
from janis_bioinformatics.tools.variantcallers import GatkSomaticVariantCaller_4_1_3
from janis_bioinformatics.tools.variantcallers.vardictsomatic_variants import (
    VardictSomaticVariantCaller,
)
from janis_pipelines.alignment.alignment import BwaAlignment
from janis_pipelines.wgs_germline.wgsgermline import WGSGermlineMultiCallers
from janis_pipelines.reference import WGS_INPUTS

INPUT_DOCS = {
    **WGS_INPUTS,
    "normal_inputs": {
        "doc": "An array of NORMAL FastqGz pairs. These are aligned separately and merged "
        "to create higher depth coverages from multiple sets of reads",
        "quality": InputQualityType.user,
        "example": [
            ["normal_R1.fastq.gz", "normal_R2.fastq.gz"],
            ["normal_R1-TOPUP.fastq.gz", "normal_R2-TOPUP.fastq.gz"],
        ],
    },
    "tumor_inputs": {
        "doc": "An array of TUMOR FastqGz pairs. These are aligned separately and merged "
        "to create higher depth coverages from multiple sets of reads",
        "quality": InputQualityType.user,
        "example": [
            ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"],
            ["tumor_R1-TOPUP.fastq.gz", "tumor_R2-TOPUP.fastq.gz"],
        ],
    },
    "normal_name": {
        "doc": "Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem",
        "quality": InputQualityType.user,
        "example": "NA12878_normal",
    },
    "tumor_name": {
        "doc": "Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem",
        "quality": InputQualityType.user,
        "example": "NA12878_tumor",
    },
    "normal_bam": {
        "doc": "Indexed NORMAL bam to call somatic variants against",
        "quality": InputQualityType.user,
        "example": "NA12878-normal.bam",
    },
    "tumor_bam": {
        "doc": "Indexed TUMOR bam to call somatic variants against",
        "quality": InputQualityType.user,
        "example": "NA12878-normal.bam",
    },
}


class WGSSomaticMultiCallers(WGSGermlineMultiCallers):
    def id(self):
        return "WGSSomaticMultiCallers"

    def friendly_name(self):
        return "WGS Somatic (Multi callers)"

    def version(self):
        return "1.4.0"

    def constructor(self):
        self.add_inputs()
        self.add_alignment_normal()
        self.add_alignment_tumor()
        self.add_bam_qc(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_gridss(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_facets(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_gatk_variantcaller(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_vardict_variantcaller(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_strelka_variantcaller(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )
        self.add_circos_plot()
        self.add_combine_variants(
            normal_bam_source=self.alignment_normal.out_bam,
            tumor_bam_source=self.alignment_tumor.out_bam,
        )

    def add_inputs(self):
        # INPUTS
        self.input("normal_inputs", Array(FastqGzPair), doc=INPUT_DOCS["normal_inputs"])
        self.input("tumor_inputs", Array(FastqGzPair), doc=INPUT_DOCS["tumor_inputs"])
        self.input("normal_name", String(), doc=INPUT_DOCS["normal_name"])
        self.input("tumor_name", String(), doc=INPUT_DOCS["tumor_name"])
        self.add_inputs_for_reference()
        self.add_inputs_for_adapter_trimming()
        self.add_inputs_for_intervals()
        self.add_inputs_for_configuration()
        self.add_inputs_facets()

    def add_inputs_for_configuration(self):
        super().add_inputs_for_configuration()
        # GATK
        self.input("gnomad", VcfTabix(), doc=INPUT_DOCS["gnomad"])
        self.input(
            "panel_of_normals",
            VcfTabix(optional=True),
            doc=INPUT_DOCS["panel_of_normals"],
        )

    def add_inputs_facets(self):
        # FACETS
        self.input("pseudo_snps", Int(optional=True))
        self.input("max_depth", Int(optional=True))
        self.input("everything", Boolean(optional=True))
        self.input("genome", String)  # Shared with CIRCOS PLOT
        self.input("cval", Int(optional=True))
        self.input("purity_cval", Int(optional=True))
        self.input("normal_depth", Int(optional=True))

    def add_alignment_normal(self):
        self.step(
            "alignment_normal",
            BwaAlignment(
                sample_name=self.normal_name,
                fastqs=self.normal_inputs,
                reference=self.reference,
                snps_dbsnp=self.snps_dbsnp,
                snps_1000gp=self.snps_1000gp,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                adapter_file=self.adapter_file,
                contaminant_file=self.contaminant_file,
            ),
        )
        self.output(
            "out_normal_R1_fastqc_reports",
            source=self.alignment_normal.out_fastqc_R1_reports,
            output_folder="reports",
        )
        self.output(
            "out_normal_R2_fastqc_reports",
            source=self.alignment_normal.out_fastqc_R2_reports,
            output_folder="reports",
        )
        self.output(
            "out_normal_bam",
            source=self.alignment_normal.out_bam,
            output_folder="bams",
            output_name=self.normal_name,
        )

    def add_alignment_tumor(self):
        self.step(
            "alignment_tumor",
            BwaAlignment(
                sample_name=self.tumor_name,
                fastqs=self.tumor_inputs,
                reference=self.reference,
                snps_dbsnp=self.snps_dbsnp,
                snps_1000gp=self.snps_1000gp,
                known_indels=self.known_indels,
                mills_indels=self.mills_indels,
                adapter_file=self.adapter_file,
                contaminant_file=self.contaminant_file,
            ),
        )
        self.output(
            "out_tumor_R1_fastqc_reports",
            source=self.alignment_tumor.out_fastqc_R1_reports,
            output_folder="reports",
        )
        self.output(
            "out_tumor_R2_fastqc_reports",
            source=self.alignment_tumor.out_fastqc_R2_reports,
            output_folder="reports",
        )
        self.output(
            "out_tumor_bam",
            source=self.alignment_tumor.out_bam,
            output_folder="bams",
            output_name=self.tumor_name,
        )

    def add_bam_qc(self, normal_bam_source, tumor_bam_source):
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
            "performance_summary_normal",
            PerformanceSummaryGenome_0_1_0(
                bam=normal_bam_source,
                genome_file=self.calculate_performancesummary_genomefile.out,
                sample_name=self.normal_name,
            ),
        )
        self.step(
            "performance_summary_tumor",
            PerformanceSummaryGenome_0_1_0(
                bam=tumor_bam_source,
                genome_file=self.calculate_performancesummary_genomefile.out,
                sample_name=self.tumor_name,
            ),
        )
        self.output(
            "out_normal_performance_summary",
            source=self.performance_summary_normal.performanceSummaryOut,
            output_folder=["performance_summary", self.normal_name],
            doc="A text file of performance summary of bam",
        )
        self.output(
            "out_tumor_performance_summary",
            source=self.performance_summary_tumor.performanceSummaryOut,
            output_folder=["performance_summary", self.tumor_name],
            doc="A text file of performance summary of bam",
        )

    def add_gridss(self, normal_bam_source, tumor_bam_source):
        # GRIDSS
        self.step(
            "vc_gridss",
            Gridss_2_6_2(
                bams=[normal_bam_source, tumor_bam_source],
                reference=self.reference,
                blacklist=self.gridss_blacklist,
            ),
        )
        self.output(
            "out_gridss_assembly",
            source=self.vc_gridss.assembly,
            output_folder=["sv", "gridss"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gridss",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Assembly returned by GRIDSS",
        )
        self.output(
            "out_variants_gridss",
            source=self.vc_gridss.out,
            output_folder=["sv", "gridss"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gridss",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Variants from the GRIDSS variant caller",
        )

    def add_facets(self, normal_bam_source, tumor_bam_source):
        # FACETS
        self.step(
            "vc_facets",
            FacetsWorkflow(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
                snps_dbsnp=self.snps_dbsnp,
                pseudo_snps=self.pseudo_snps,
                max_depth=self.max_depth,
                everything=self.everything,
                genome=self.genome,
                cval=self.cval,
                purity_cval=self.purity_cval,
                normal_depth=self.normal_depth,
            ),
        )
        self.output(
            "out_facets_summary",
            source=self.vc_facets.out_summary,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        )
        self.output(
            "out_facets_purity_png",
            source=self.vc_facets.out_purity_png,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_purity",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_purity_seg",
            source=self.vc_facets.out_purity_seg,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_purity",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_purity_rds",
            source=self.vc_facets.out_purity_rds,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_purity",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_png",
            source=self.vc_facets.out_hisens_png,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_seg",
            source=self.vc_facets.out_hisens_seg,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_rds",
            source=self.vc_facets.out_hisens_rds,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_arm_level",
            source=self.vc_facets.out_arm_level,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}.arm_level",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_gene_level",
            source=self.vc_facets.out_gene_level,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}.gene_level",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_qc",
            source=self.vc_facets.out_qc,
            output_folder=["cnv", "facets"],
            output_name=StringFormatter(
                "{tumour}--{normal}.qc",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),

    def add_gatk_variantcaller(self, normal_bam_source, tumor_bam_source):
        # GATK
        if "generate_gatk_intervals" in self.step_nodes:
            generated_intervals = self.generate_gatk_intervals.out_regions
        else:
            generated_intervals = self.step(
                "generate_gatk_intervals",
                GenerateIntervalsByChromosome(reference=self.reference),
                when=self.gatk_intervals.is_null(),
            ).out_regions
        intervals = FirstOperator([self.gatk_intervals, generated_intervals])
        recal_ins = {
            "reference": self.reference,
            "intervals": intervals,
            "snps_dbsnp": self.snps_dbsnp,
            "snps_1000gp": self.snps_1000gp,
            "known_indels": self.known_indels,
            "mills_indels": self.mills_indels,
        }
        self.step(
            "bqsr_normal",
            GATKBaseRecalBQSRWorkflow_4_1_3(bam=normal_bam_source, **recal_ins),
            scatter="intervals",
        )
        self.step(
            "bqsr_tumor",
            GATKBaseRecalBQSRWorkflow_4_1_3(bam=tumor_bam_source, **recal_ins),
            scatter="intervals",
        )
        self.step(
            "vc_gatk",
            GatkSomaticVariantCaller_4_1_3(
                normal_bam=self.bqsr_normal.out,
                tumor_bam=self.bqsr_tumor.out,
                normal_name=self.normal_name,
                intervals=intervals,
                reference=self.reference,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
            ),
            scatter=["intervals", "normal_bam", "tumor_bam"],
        )
        self.step(
            "vc_gatk_merge",
            BcfToolsConcat_1_9(vcf=self.vc_gatk.out.as_type(Array(Vcf))),
        )
        self.step(
            "vc_gatk_sort_combined",
            BcfToolsSort_1_9(vcf=self.vc_gatk_merge.out.as_type(CompressedVcf)),
        )
        self.step(
            "vc_gatk_uncompressvcf",
            UncompressArchive(
                file=self.vc_gatk_sort_combined.out,
                force=True
            ),
        )
        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_sort_combined.out,
            output_folder="variants",
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gatk",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Merged variants from the GATK caller",
        )
        self.output(
            "out_variants_gatk_split",
            source=self.vc_gatk.out,
            output_folder=["variants", "GatkByInterval"],
            doc="Unmerged variants from the GATK caller (by interval)",
        )

    def add_strelka_variantcaller(self, normal_bam_source, tumor_bam_source):
        self.step("generate_manta_config", GenerateMantaConfig())
        self.step(
            "vc_strelka",
            IlluminaSomaticVariantCaller(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                intervals=self.strelka_intervals,
                reference=self.reference,
                manta_config=self.generate_manta_config.out,
            ),
        )
        self.step("vc_strelka_compress", BGZipLatest(file=self.vc_strelka.out))
        self.output(
            "out_variants_strelka",
            source=self.vc_strelka_compress.out.as_type(CompressedVcf),
            output_folder=["variants"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_strelka",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Variants from the Strelka variant caller",
        )
        self.output(
            "out_variants_manta_somatic",
            source=self.vc_strelka.tumor_sv,
            output_folder=["sv", "manta"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_manta",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="SV variants from the Manta caller",
        )

    def add_vardict_variantcaller(self, normal_bam_source, tumor_bam_source):
        self.step(
            "generate_vardict_headerlines",
            GenerateVardictHeaderLines(reference=self.reference),
        )
        self.step(
            "vc_vardict",
            VardictSomaticVariantCaller(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                normal_name=self.normal_name,
                tumor_name=self.tumor_name,
                header_lines=self.generate_vardict_headerlines.out,
                intervals=self.vardict_intervals,
                reference=self.reference,
                allele_freq_threshold=self.allele_freq_threshold,
                minMappingQual=self.minMappingQual,
                filter=self.filter,
            ),
            scatter="intervals",
        )
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
            output_folder=["variants"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_vardict",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Merged variants from the VarDict caller",
        )
        self.output(
            "out_variants_vardict_split",
            source=self.vc_vardict.out,
            output_folder=["variants", "VardictByInterval"],
            doc="Unmerged variants from the GATK caller (by interval)",
        )

    def add_circos_plot(self):
        self.step(
            "circos_plot",
            CircosPlot_0_1_2(
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
                facets_file=self.vc_facets.out_hisens_rds,
                sv_file=self.vc_strelka.tumor_sv.assert_not_null(),
                genome=self.genome,
            ),
        )
        self.output(
            "out_circos_plot",
            source=self.circos_plot.out,
            output_folder="circos_plot",
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_circos_plot",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
        )

    def add_combine_variants(self, normal_bam_source, tumor_bam_source):
        self.step(
            "combine_variants",
            CombineVariants_0_0_8(
                normal=self.normal_name,
                tumor=self.tumor_name,
                vcfs=[
                    self.vc_gatk_uncompressvcf.out.as_type(Vcf),
                    self.vc_strelka.out,
                    self.vc_vardict_uncompress_for_combine.out.as_type(Vcf),
                ],
                type="somatic",
                columns=["AD", "DP", "GT"],
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
            AddBamStatsSomatic_0_1_0(
                normal_id=self.normal_name,
                tumor_id=self.tumor_name,
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                vcf=self.combined_uncompress.out.as_type(Vcf),
                reference=self.reference,
            ),
        )

        self.output(
            "out_variants_combined_bamstats",
            source=self.combined_addbamstats.out,
            output_folder=["variants"],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_combined",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Combined variants from GATK, VarDict and Strelka callers",
        )

    def bind_metadata(self):
        meta: WorkflowMetadata = super().bind_metadata() or self.metadata

        meta.keywords = [
            "wgs",
            "cancer",
            "somatic",
            "variants",
            "gatk",
            "vardict",
            "strelka",
            "gridss",
            "facets",
        ]
        meta.contributors = ["Michael Franklin", "Richard Lupat", "Jiaan Yu"]
        meta.dateCreated = date(2018, 12, 24)
        meta.dateUpdated = date(2022, 3, 1)
        meta.short_documentation = "A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2."
        meta.documentation = """\
This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs:

- Takes raw sequence data in the FASTQ format;
- align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate somatic variant callers (GATK / Strelka / VarDict);
- Outputs the final variants in the VCF format.

**Resources**

This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:

- https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

This pipeline expects the assembly references to be as they appear in that storage \
    (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.
"""

    def tests(self) -> Optional[List[TTestCase]]:
        parent_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics"
        brca1_test_data = f"{parent_dir}/brca1_test/test_data"

        # To do: Vcf.basic_test to remove headers and md5
        return [
            TTestCase(
                name="brca1",
                input={
                    "normal_inputs": [
                        [
                            f"{brca1_test_data}/NA24385-BRCA1_R1.fastq.gz",
                            f"{brca1_test_data}/NA24385-BRCA1_R2.fastq.gz",
                        ]
                    ],
                    "normal_name": "NA24385-BRCA1",
                    "tumor_inputs": [
                        [
                            f"{brca1_test_data}/NA12878-NA24385-mixture-BRCA1_R1.fastq.gz",
                            f"{brca1_test_data}/NA12878-NA24385-mixture-BRCA1_R2.fastq.gz",
                        ]
                    ],
                    "tumor_name": "NA12878-NA24385-mixture-BRCA1",
                    "reference": f"{brca1_test_data}/Homo_sapiens_assembly38.chr17.fasta",
                    "gridss_blacklist": f"{brca1_test_data}/consensusBlacklist.hg38.chr17.bed",
                    "gnomad": f"{brca1_test_data}/af-only-gnomad.hg38.BRCA1.vcf.gz",
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
                    "genome": f"hg38",
                    "vc_strelka_manta_runtime_cpu": 1,
                    "vc_strelka_strelka_runtime_cpu": 1,
                    "vc_strelka_manta_runtime_memory": 5,
                    "vc_strelka_strelka_runtime_memory": 3,
                },
                output=Array.array_wrapper(
                    [ZipFile.basic_test("out_normal_R1_fastqc_reports", 430000)]
                )
                + Array.array_wrapper(
                    [ZipFile.basic_test("out_tumor_R1_fastqc_reports", 430000)]
                )
                + Array.array_wrapper(
                    [ZipFile.basic_test("out_normal_R2_fastqc_reports", 430000)]
                )
                + Array.array_wrapper(
                    [ZipFile.basic_test("out_tumor_R2_fastqc_reports", 430000)]
                )
                + TextFile.basic_test(
                    "out_normal_performance_summary",
                    950,
                    md5="e3205735e5fe8c900f05050f8ed73f19",
                )
                + TextFile.basic_test(
                    "out_tumor_performance_summary",
                    950,
                    md5="122bfa2ece90c0f030015feba4ba7d84",
                )
                + BamBai.basic_test("out_normal_bam", 3260000, 49000)
                + BamBai.basic_test("out_tumor_bam", 3340000, 49000)
                + Bam.basic_test("out_gridss_assembly", 60000)
                + Vcf.basic_test("out_variants_gridss", 90000)
                + File.basic_test(
                    "out_facets_summary", 450, "c9b543b6e62003b8b6409db1a4e78248"
                )
                + File.basic_test("out_facets_purity_png", 40000)
                + File.basic_test("out_facets_purity_seg", 100)
                + File.basic_test("out_facets_purity_rds", 7000)
                + File.basic_test("out_facets_hisens_png", 40000)
                + File.basic_test("out_facets_hisens_seg", 100)
                + File.basic_test("out_facets_hisens_rds", 7000)
                + CompressedVcf.basic_test("out_variants_gatk", 9000, 149)
                + Array.array_wrapper(
                    [Vcf.basic_test("out_variants_gatk_split", 34000, 147)]
                )
                + CompressedVcf.basic_test("out_variants_vardict", 13000, 189)
                + Array.array_wrapper(
                    [Vcf.basic_test("out_variants_vardict_split", 55000, 187)]
                )
                + CompressedVcf.basic_test("out_variants_strelka", 7000)
                + VcfTabix.basic_test("out_variants_manta_somatic", 1400, 70, 35)
                + File.basic_test("out_circos_plot", 60000)
                + Vcf.basic_test("out_variants_combined_bamstats", 91000),
            )
        ]


if __name__ == "__main__":
    # import os.path

    # w = WGSSomaticMultiCallers()
    # args = {
    #     "to_console": False,
    #     "to_disk": True,
    #     "validate": True,
    #     "export_path": os.path.join(
    #         os.path.dirname(os.path.realpath(__file__)), "{language}"
    #     ),
    # }
    # w.translate("cwl", **args)
    # w.translate("wdl", **args)
    WGSSomaticMultiCallers().translate("wdl")
