
from datetime import datetime

from janis_core.tool.test_classes import TTestCase
from janis_core import ToolMetadata

from ..types import BamBai, FastaWithDict, VcfTabix, Bed
from ..tools import Gatk4BaseRecalibrator_4_1_3, Gatk4ApplyBqsr_4_1_3
from ..workflows import BioinformaticsWorkflow


class GATKBaseRecalBQSRWorkflow_4_1_3(BioinformaticsWorkflow):
    def id(self):
        return "GATKBaseRecalBQSRWorkflow"

    def friendly_name(self):
        return "GATK Base Recalibration on Bam"

    def version(self):
        return "4.1.3"

    def tool_provider(self):
        return "common"

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
        self.input("snps_1000gp", VcfTabix)
        self.input("known_indels", VcfTabix)
        self.input("mills_indels", VcfTabix)

        self.step(
            "base_recalibrator",
            Gatk4BaseRecalibrator_4_1_3(
                bam=self.bam,
                intervals=self.intervals,
                reference=self.reference,
                knownSites=[
                    self.snps_dbsnp,
                    self.snps_1000gp,
                    self.known_indels,
                    self.mills_indels,
                ],
            ),
        )
        self.step(
            "apply_bqsr",
            Gatk4ApplyBqsr_4_1_3(
                bam=self.bam,
                intervals=self.intervals,
                recalFile=self.base_recalibrator.out,
                reference=self.reference,
            ),
        )
        self.output("out", source=self.apply_bqsr.out)

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Jiaan Yu"],
            dateCreated=datetime(2020, 6, 12),
            dateUpdated=datetime(2020, 6, 12),
            documentation="",
        )

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "snps_dbsnp": f"{remote_dir}/Homo_sapiens_assembly38.dbsnp138.BRCA1.vcf.gz",
                    "snps_1000gp": f"{remote_dir}/1000G_phase1.snps.high_confidence.hg38.BRCA1.vcf.gz",
                    "known_indels": f"{remote_dir}/Homo_sapiens_assembly38.known_indels.BRCA1.vcf.gz",
                    "mills_indels": f"{remote_dir}/Mills_and_1000G_gold_standard.indels.hg38.BRCA1.vcf.gz",
                    "intervals": f"{remote_dir}/BRCA1.hg38.bed",
                },
                output=BamBai.basic_test(
                    "out",
                    2600000,
                    21000,
                    f"{remote_dir}/NA12878-BRCA1.recalibrated.flagstat",
                ),
            )
        ]
