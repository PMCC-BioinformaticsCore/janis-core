

from datetime import datetime
from janis_core import Array, Boolean, Int, String, ToolMetadata
from janis_core.tool.test_classes import TTestCase

from ..types import BamBai
from ..tools import Gatk4MarkDuplicates_4_1_3
from ..tools import Gatk4MergeSamFiles_4_1_3
from . import BioinformaticsWorkflow


class MergeAndMarkBams_4_1_3(BioinformaticsWorkflow):
    def id(self):
        return "mergeAndMarkBams"

    def friendly_name(self):
        return "Merge and Mark Duplicates"

    def version(self):
        return "4.1.3"

    def tool_provider(self):
        return "common"

    def constructor(self):

        self.input("bams", Array(BamBai()))
        self.input("createIndex", Boolean, default=True)
        self.input("maxRecordsInRam", Int, default=5000000)
        self.input("sampleName", String(optional=True))

        self.step(
            "mergeSamFiles",
            Gatk4MergeSamFiles_4_1_3(
                bams=self.bams,
                useThreading=True,
                createIndex=self.createIndex,
                maxRecordsInRam=self.maxRecordsInRam,
                validationStringency="SILENT",
                sampleName=self.sampleName,
            ),
        )

        self.step(
            "markDuplicates",
            Gatk4MarkDuplicates_4_1_3(
                bam=self.mergeSamFiles.out,
                createIndex=self.createIndex,
                maxRecordsInRam=self.maxRecordsInRam,
                outputPrefix=self.sampleName,
            ),
        )
        self.output("out", source=self.markDuplicates.out)

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 2, 19),
            dateUpdated=datetime(2020, 11, 6),
            documentation="",
        )

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "bams": [
                        f"{remote_dir}/NA12878-BRCA1.sorted.bam",
                    ],
                    "maxRecordsInRam": 5000000,
                    "createIndex": True,
                    "mergeSamFiles_useThreading": True,
                    "mergeSamFiles_validationStringency": "SILENT",
                },
                output=BamBai.basic_test(
                    "out",
                    2829000,
                    3780,
                   f"{remote_dir}/NA12878-BRCA1.markduped.bam.flagstat",
                ),
            ),
            TTestCase(
                name="minimal",
                input={
                    "bams": [
                        f"{remote_dir}/NA12878-BRCA1.sorted.bam",
                    ],
                    "maxRecordsInRam": 5000000,
                    "createIndex": True,
                    "mergeSamFiles_useThreading": True,
                    "mergeSamFiles_validationStringency": "SILENT",
                },
                output=self.minimal_test(),
            ),
        ]


if __name__ == "__main__":
    MergeAndMarkBams_4_1_3().translate("wdl")
