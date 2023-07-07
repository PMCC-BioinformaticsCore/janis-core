

import os
from abc import ABC
import datetime

from janis_core import ToolInput, ToolOutput, File, Filename, InputSelector, String, CommandTool
from janis_core.tool.test_classes import TTestCase
from .bioinformaticstool import BioinformaticsTool
from ..types import Vcf
from .petermacutils import PeterMacUtils_0_0_7




class AddBamStatsBase(BioinformaticsTool, ABC):
    def tool(self):
        return "addBamStats"

    def friendly_name(self):
        return "Add Bam Statistics to Vcf"

    def base_command(self):
        return "add_bam_stats.py"

    def inputs(self):
        return [
            *self.additional_inputs,
            ToolInput("inputVcf", Vcf(), prefix="-i", doc="input vcf"),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf", suffix=".addbamstats"),
                prefix="-o",
                doc="output vcf name",
            ),
            ToolInput(
                "type",
                String(),
                prefix="--type",
                doc="must be either germline or somatic",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Vcf(), glob=InputSelector("outputFilename"))]

    def bind_metadata(self):
        self.metadata.dateCreated = datetime.datetime(2020, 5, 20)
        self.metadata.dateUpdated = datetime.datetime(2020, 5, 20)
        self.metadata.contributors = ["Jiaan Yu"]
        self.metadata.documentationUrl = (
            "https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils"
        )
        self.metadata.documentation = """usage: add_bam_stats.py [-h] -i I -o O --type {germline,somatic}
                        [--mpileup MPILEUP] [--normal_mpileup NORMAL_MPILEUP]
                        [--tumor_mpileup TUMOR_MPILEUP]
                        [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

Get stats from bam file and write to vcf

required arguments:
  -i I                  input vcf
  -o O                  output vcf
  --type {germline,somatic}
                        must be either germline or somatic
  --mpileup MPILEUP     mpileup file extracted from bam file
  --normal_mpileup NORMAL_MPILEUP
                        mpileup file extracted from the normal sample bam,
                        required if input is somatic vcf
  --tumor_mpileup TUMOR_MPILEUP
                        mpileup file extracted from the tumor sample, required
                        if input is somatic vcf
  --normal_id NORMAL_ID
                        Normal sample id, required if input is somatic vcf
  --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

optional arguments:
  -h, --help            show this help message and exit
        """

    additional_inputs = [
        ToolInput(
            "mpileup",
            File(optional=True),
            prefix="--mpileup",
            doc="mpileup file extracted from bam file",
        ),
        ToolInput(
            "normalMpileup",
            File(optional=True),
            prefix="--normal_mpileup",
            doc="mpileup file extracted from the normal sample bam, required if input is somatic vcf",
        ),
        ToolInput(
            "tumorMpileup",
            File(optional=True),
            prefix="--tumor_mpileup",
            doc="mpileup file extracted from the tumor sample bam, required if input is somatic vcf",
        ),
        ToolInput(
            "normalID",
            String(optional=True),
            prefix="--normal_id",
            doc="normal sample id, required if input is somatic vcf",
        ),
        ToolInput(
            "tumorID",
            String(optional=True),
            prefix="--tumor_id",
            doc="tumor sample id, required if input is somatic vcf",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "inputVcf": f"{remote_dir}/NA12878-BRCA1.sorted.uncompressed.stdout",
                    "mpileup": f"{remote_dir}/NA12878-BRCA1.mpileup.stdout",
                    "type": "germline",
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



# NOTE: disabled as this is for dev work only
# class AddBamStats_dev(AddBamStatsBase, PeterMacUtils_dev):
#     pass


class AddBamStats_0_0_7(AddBamStatsBase, PeterMacUtils_0_0_7):
    pass


AddBamStatsLatest = AddBamStats_0_0_7
