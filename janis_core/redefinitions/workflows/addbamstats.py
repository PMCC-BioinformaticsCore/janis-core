
import os
from datetime import datetime
from janis_core import WorkflowBuilder, WorkflowMetadata
from janis_core.tool.test_classes import TTestCase
from janis_core import String

from .bioinformaticsworkflow import BioinformaticsWorkflow
from ..types import Vcf, BamBai, FastaWithDict
from ..tools import AddBamStatsLatest
from ..tools import SamToolsMpileupLatest


class AddBamStatsGermline_0_1_0(BioinformaticsWorkflow):
    def id(self) -> str:
        return "AddBamStatsGermline"

    def friendly_name(self):
        return "Annotate Bam Stats to Germline Vcf Workflow"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def bind_metadata(self):
        return WorkflowMetadata(
            version="v0.1.0",
            contributors=["Jiaan Yu"],
            dateCreated=datetime(2020, 6, 4),
            dateUpdated=datetime(2020, 8, 10),
        )

    def constructor(self):

        # self.input("sample_name", String)
        self.input("bam", BamBai)
        self.input("vcf", Vcf)
        self.input("reference", FastaWithDict)

        self.step(
            "samtoolsmpileup",
            SamToolsMpileupLatest(
                bam=self.bam,
                positions=self.vcf,
                reference=self.reference,
                countOrphans=True,
                noBAQ=True,
                minBQ=0,
                maxDepth=10000,
            ),
        )

        self.step(
            "addbamstats",
            AddBamStatsLatest(
                inputVcf=self.vcf,
                mpileup=self.samtoolsmpileup.out,
                type="germline",
            ),
        )

        self.output("out", source=self.addbamstats.out, output_name="addbasmtats.vcf")

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bam": f"{remote_dir}/NA12878-BRCA1.markduped.bam",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                    "vcf": f"{remote_dir}/NA12878-BRCA1.sorted.uncompressed.stdout",
                    "samtoolsmpileup_countOrphans": True,
                    "samtoolsmpileup_noBAQ": True,
                    "samtoolsmpileup_maxDepth": 10000,
                    "samtoolsmpileup_minBQ": 0,
                    "addbamstats_type": "germline",
                },
                output=Vcf.basic_test(
                    "out",
                    69225,
                    230,
                    ["GATKCommandLine"],
                    "db09c6c37c52771bd058e32d5c6b94c1",
                ),
            )
        ]
