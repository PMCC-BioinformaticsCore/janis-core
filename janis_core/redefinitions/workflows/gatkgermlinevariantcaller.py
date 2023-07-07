
from datetime import date

from janis_core.tool.test_classes import TTestCase

from ..types import FastaWithDict, BamBai, VcfTabix, Bed, Vcf
from ..tools import UncompressArchive
from ..tools import Gatk4SplitReads_4_1_3
from ..tools import Gatk4HaplotypeCaller_4_1_3
from ..tools import SplitMultiAllele
from .bioinformaticsworkflow import BioinformaticsWorkflow


class GatkGermlineVariantCaller_4_1_3(BioinformaticsWorkflow):
    def id(self):
        return "GATK4_GermlineVariantCaller"

    def friendly_name(self):
        return "GATK4 Germline Variant Caller"

    def tool_provider(self):
        return "Variant Callers"

    def bind_metadata(self):
        self.metadata.version = "4.1.3.0"
        self.metadata.dateCreated = date(2019, 9, 1)
        self.metadata.dateUpdated = date(2019, 9, 13)

        self.metadata.contributors = ["Michael Franklin", "Jiaan"]
        self.metadata.keywords = ["variants", "gatk", "gatk4", "variant caller"]
        self.metadata.documentation = """
        This is a VariantCaller based on the GATK Best Practice pipelines. It uses the GATK4 toolkit, specifically 4.1.3.

        It has the following steps:

        1. Split Bam based on intervals (bed)
        2. HaplotypeCaller
        3. SplitMultiAllele
                """.strip()

    def constructor(self):

        self.input("bam", BamBai)
        self.input(
            "intervals",
            Bed(optional=True),
            doc="This optional interval supports processing by regions. If this input resolves "
            "to null, then GATK will process the whole genome per each tool's spec",
        )
        self.input("reference", FastaWithDict)
        self.input("snps_dbsnp", VcfTabix)

        self.step(
            "split_bam",
            Gatk4SplitReads_4_1_3(outputFilename='out', bam=self.bam, intervals=self.intervals),
        )

        self.step(
            "haplotype_caller",
            Gatk4HaplotypeCaller_4_1_3(
                inputRead=self.split_bam.out,
                intervals=self.intervals,
                reference=self.reference,
                dbsnp=self.snps_dbsnp,
                pairHmmImplementation="LOGLESS_CACHING",
            ),
        )
        self.step("uncompressvcf", UncompressArchive(file=self.haplotype_caller.out, force=True))
        self.step(
            "splitnormalisevcf",
            SplitMultiAllele(
                vcf=self.uncompressvcf.out.as_type(Vcf), reference=self.reference
            ),
        )

        self.output("variants", source=self.haplotype_caller.out)
        self.output("out_bam", source=self.haplotype_caller.bam)
        self.output("out", source=self.splitnormalisevcf.out)

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.recalibrated.bam",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "snps_dbsnp": f"{remote_dir}/Homo_sapiens_assembly38.dbsnp138.BRCA1.vcf.gz",
                    "haplotype_caller_pairHmmImplementation": "LOGLESS_CACHING",
                },
                output=Vcf.basic_test(
                    "out",
                    51000,
                    221,
                    ["GATKCommandLine"],
                    "5e48624cb5ef379a7d6d39cec44bc856",
                ),
            )
        ]


if __name__ == "__main__":
    vc = GatkGermlineVariantCaller_4_1_3().translate("wdl", to_console=True)
    # print(vc.translate("cwl"))
