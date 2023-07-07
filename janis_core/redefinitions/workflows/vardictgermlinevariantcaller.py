

from datetime import datetime
from janis_core import File, String, Float, Int, WorkflowMetadata

from ..types import FastaWithDict, BamBai, Bed
from ..tools import SplitMultiAllele
from ..tools import BGZipLatest
from ..tools import VcfToolsvcftoolsLatest
from ..tools import BcfToolsAnnotate_1_5
from ..tools import TabixLatest
from ..tools import VarDictGermline_1_6_0
from ..tools import TrimIUPAC_0_0_5
from .bioinformaticsworkflow import BioinformaticsWorkflow


class VardictGermlineVariantCaller(BioinformaticsWorkflow):
    def id(self):
        return "vardictGermlineVariantCaller"

    def friendly_name(self):
        return "Vardict Germline Variant Caller"

    def tool_provider(self):
        return "Variant Callers"

    def version(self):
        return "v0.1.1"

    def constructor(self):

        self.input("bam", BamBai)
        self.input("intervals", Bed)
        self.input("sample_name", String)
        self.input("header_lines", File)
        self.input("reference", FastaWithDict)

        # vardict options
        self.input("allele_freq_threshold", Float, default=0.05)
        self.input("minMappingQual", Int(optional=True))
        self.input("filter", String(optional=True))

        self.step(
            "vardict",
            VarDictGermline_1_6_0(
                intervals=self.intervals,
                bam=self.bam,
                reference=self.reference,
                sampleName=self.sample_name,
                var2vcfSampleName=self.sample_name,
                alleleFreqThreshold=self.allele_freq_threshold,
                var2vcfAlleleFreqThreshold=self.allele_freq_threshold,
                vcfFormat=True,
                chromColumn=1,
                regStartCol=2,
                geneEndCol=3,
                threads=4,
                minMappingQual=self.minMappingQual,
                filter=self.filter,
            ),
        )
        self.step(
            "annotate",
            BcfToolsAnnotate_1_5(vcf=self.vardict.out, headerLines=self.header_lines),
        )
        self.step("compressvcf", BGZipLatest(file=self.annotate.out, stdout=True))
        self.step("tabixvcf", TabixLatest(inp=self.compressvcf.out))

        self.step(
            "splitnormalisevcf",
            SplitMultiAllele(vcf=self.annotate.out, reference=self.reference),
        )
        self.step("trim", TrimIUPAC_0_0_5(vcf=self.splitnormalisevcf.out))
        self.step(
            "filterpass",
            VcfToolsvcftoolsLatest(
                vcf=self.trim.out,
                removeFileteredAll=True,
                recode=True,
                recodeINFOAll=True,
            ),
        )

        self.output("variants", source=self.tabixvcf.out)
        self.output("out", source=self.filterpass.out)

    def bind_metadata(self):
        return WorkflowMetadata(
            contributors=["Michael Franklin", "Jiaan Yu"],
            dateCreated=datetime(2019, 3, 28),
            dateUpdated=datetime(2021, 3, 5),
            documentation="",
        )


if __name__ == "__main__":
    v = VardictGermlineVariantCaller()
    v.translate("wdl", with_resource_overrides=False)
    # print(v.generate_resources_file("wdl", { "CaptureType": "targeted" }))
