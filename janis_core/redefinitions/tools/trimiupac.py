
from datetime import datetime
from abc import ABC
from typing import Dict, Any, List

from janis_core import (
    CaptureType,
    ToolInput,
    Filename,
    ToolOutput,
    InputSelector,
    ToolMetadata,
)
from janis_core import get_value_for_hints_and_ordered_resource_tuple

from .petermacutils import PeterMacUtils_0_0_4
from .petermacutils import PeterMacUtils_0_0_5

from ..types import Vcf
from .bioinformaticstool import BioinformaticsTool


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
            CaptureType.CHROMOSOME: 2,
            CaptureType.EXOME: 2,
            CaptureType.THIRTYX: 4,
            CaptureType.NINETYX: 4,
            CaptureType.THREEHUNDREDX: 4,
        },
    )
]


class TrimIUPACBase(BioinformaticsTool, ABC):
    def tool(self) -> str:
        return "trimIUPAC"

    def friendly_name(self) -> str:
        return "Trim IUPAC Bases"

    def base_command(self):
        return "trimIUPAC.py"

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 1

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 1

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput(
                "vcf", Vcf(), position=0, doc="The VCF to remove the IUPAC bases from"
            ),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf", suffix=".trimmed"),
                position=2,
            ),
        ]

    def outputs(self) -> List[ToolOutput]:
        return [ToolOutput("out", Vcf(), InputSelector("outputFilename"))]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 5, 30),
            dateUpdated=datetime(2019, 12, 8),
            documentation="",
        )


class TrimIUPAC_0_0_4(TrimIUPACBase, PeterMacUtils_0_0_4):
    pass


class TrimIUPAC_0_0_5(TrimIUPACBase, PeterMacUtils_0_0_5):
    pass


TrimIUPACLatest = TrimIUPAC_0_0_5
