

from datetime import datetime
from typing import List, Dict, Any
from janis_core import (
    get_value_for_hints_and_ordered_resource_tuple,
    ToolMetadata,
    UnionType,
)
from janis_core.tool.test_classes import TTestCase

from ..types import FastaWithDict, CompressedVcf, Vcf
from .bioinformaticstool import BioinformaticsTool
from janis_core import (
    ToolOutput,
    ToolInput,
    Filename,
    ToolArgument,
    InputSelector,
    CaptureType,
)


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
            CaptureType.CHROMOSOME: 8,
            CaptureType.EXOME: 8,
            CaptureType.THIRTYX: 8,
            CaptureType.NINETYX: 12,
            CaptureType.THREEHUNDREDX: 16,
        },
    )
]


class SplitMultiAllele(BioinformaticsTool):
    def tool(self):
        return "SplitMultiAllele"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def friendly_name(self):
        return "Split Multiple Alleles"

    def base_command(self):
        return None

    def container(self):
        return "heuermh/vt"  # "SEE (mfranklin's notes in) DOCUMENTATION"

    def version(self):
        return "v0.5772"

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

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput(
                "vcf", UnionType(Vcf, CompressedVcf), position=1, shell_quote=False
            ),
            ToolInput(
                "reference", FastaWithDict(), prefix="-r", position=4, shell_quote=False
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("vcf", remove_file_extension=True),
                    extension=".vcf",
                    suffix=".norm",
                ),
                position=6,
                prefix="-o",
                shell_quote=False,
            ),
        ]

    def arguments(self):
        return [
            ToolArgument("vt decompose -s ", position=0, shell_quote=False),
            ToolArgument("| vt normalize -n -q - ", position=2, shell_quote=False),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", Vcf(), glob=InputSelector("outputFilename"))]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin", "Jiaan Yu"],
            dateCreated=datetime(2019, 1, 18),
            dateUpdated=datetime(2020, 11, 6),
            documentation="",
        )

    def doc(self):
        return """
    Use vt to split multiallelic variants, and left-most align normalisation.
    Original command:
    vt decompose -s $input.vcf | vt normalize -n -q - -r $reference -o $output.vcf

    ========
    VT decompose documentation:
        options : -s  smart decomposition [false]
              -d  debug [false]
              -f  filter expression []
              -o  output VCF file [-]
              -I  file containing list of intervals []
              -i  intervals []
              -?  displays help
              
    VT normalize documentation:
        options : -o  output VCF file [-]
              -d  debug [false]
              -q  do not print options and summary [false]
              -m  warns but does not exit when REF is inconsistent
                  with masked reference sequence for non SNPs.
                  This overides the -n option [false]
              -n  warns but does not exit when REF is inconsistent
                  with reference sequence for non SNPs [false]
              -f  filter expression []
              -w  window size for local sorting of variants [10000]
              -I  file containing list of intervals []
              -i  intervals []
              -r  reference sequence fasta file []
              -?  displays help 
        """.strip()

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "vcf": f"{remote_dir}/NA12878-BRCA1.haplotype_uncompressed.stdout",
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                },
                output=Vcf.basic_test(
                    "out",
                    51462,
                    221,
                    ["GATKCommandLine"],
                    "5e48624cb5ef379a7d6d39cec44bc856",
                ),
            )
        ]


if __name__ == "__main__":
    print(SplitMultiAllele().help())
