

from abc import ABC
import datetime

from janis_core.tool.test_classes import TTestCase
from janis_core import ToolInput, ToolOutput, File, Filename, InputSelector, Boolean

from ..types import Csv, TextFile
from .bioinformaticstool import BioinformaticsTool
from .petermacutils import PeterMacUtils_0_0_7


class PerformanceSummaryBase(BioinformaticsTool, ABC):
    def tool(self):
        return "performanceSummary"

    def friendly_name(self):
        return "Performance Summary"

    def base_command(self):
        return "performance_summary.py"

    def inputs(self):
        return [
            ToolInput(
                "flagstat",
                File(),
                prefix="--flagstat",
                doc="output of samtools flagstat on bam",
            ),
            ToolInput(
                "collectInsertSizeMetrics",
                File,
                prefix="--collect_insert_metrics",
                doc="output of CollectInsertMetrics (GATK or Picard) on bam",
            ),
            ToolInput(
                "coverage",
                File(),
                prefix="--coverage",
                doc="output of bedtools coverageBed for targeted bam; bedtools genomeCoverageBed for whole genome bam",
            ),
            ToolInput(
                "outputPrefix",
                Filename(),
                prefix="-o",
                doc="prefix of output summary csv",
            ),
            *self.additional_args,
        ]

    def outputs(self):
        return [ToolOutput("out", Csv(), glob=InputSelector("outputPrefix") + ".csv")]

    additional_args = [
        ToolInput(
            "targetFlagstat",
            File(optional=True),
            prefix="--target_flagstat",
            doc="output of samtools flagstat of bam target on target bed. Only specified for targeted bam",
        ),
        ToolInput(
            "rmdupFlagstat",
            File(optional=True),
            prefix="--rmdup_flagstat",
            doc="output of samtools flagstat of removed duplicates bam. File to be used to extract mapping infomation if specified, instead of the --flagstat file.",
        ),
        ToolInput(
            "genome",
            Boolean(optional=True),
            prefix="--genome",
            doc="calculate statistics for whole genome data.--target_flagstat must not be speicified",
        ),
    ]

    def bind_metadata(self):
        self.metadata.dateCreated = datetime.datetime(2020, 4, 3)
        self.metadata.dateUpdated = datetime.datetime(2020, 4, 3)
        self.metadata.contributors = ["Jiaan Yu"]
        self.metadata.documentation = """usage: performance_summary.py [-h] --flagstat FLAGSTAT
                              --collect_insert_metrics COLLECT_INSERT_METRICS
                              --coverage COVERAGE -o O
                              [--target_flagstat TARGET_FLAGSTAT]
                              [--rmdup_flagstat RMDUP_FLAGSTAT] [--genome]

Performance summary of bam

required arguments:
  --flagstat FLAGSTAT   output of samtools flagstat on bam
  --collect_insert_metrics COLLECT_INSERT_METRICS
                        output of CollectInsertMetrics (GATK or Picard) on bam
  --coverage COVERAGE   output of bedtools coverageBed for targeted bam;
                        bedtools genomeCoverageBed for whole genome bam
  -o O                  output summary csv name

optional arguments:
  -h, --help            show this help message and exit
  --target_flagstat TARGET_FLAGSTAT
                        output of samtools flagstat of bam target on target
                        bed. Only specified for targeted bam
  --rmdup_flagstat RMDUP_FLAGSTAT
                        output of samtools flagstat of removed duplicates bam.
                        File to be used to extract mapping infomation if
                        specified, instead of the --flagstat file.
  --genome              calculate statistics for whole genome data.
                        --target_flagstat must not be speicified
        """
        self.metadata.documentationUrl = (
            "https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance"
        )

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "flagstat": f"{remote_dir}/NA12878-BRCA1.markduped.bam.flagstat",
                    "collectInsertSizeMetrics": f"{remote_dir}/NA12878-BRCA1.markduped.metrics.txt",
                    "coverage": f"{remote_dir}/NA12878-BRCA1.genomeCoverageBed.stdout",
                    "rmdupFlagstat": f"{remote_dir}/NA12878-BRCA1.markduped.bam.bam.flagstat",
                    "genome": True,
                },
                output=TextFile.basic_test(
                    tag="out",
                    min_size=948,
                    line_count=2,
                    md5="575354942cfb8d0367725f9020181443",
                    expected_file_path=f"{remote_dir}/NA12878-BRCA1_performance_summary.csv",
                ),
            )
        ]


class PerformanceSummary_0_0_7(PerformanceSummaryBase, PeterMacUtils_0_0_7):
    pass


PerformanceSummaryLatest = PerformanceSummary_0_0_7
# PerformanceSummaryLatest = PerformanceSummary_dev

