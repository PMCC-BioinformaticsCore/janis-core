

from typing import Optional
from janis_bioinformatics.data_types import BamBai, FileTabix
from janis_unix import Gunzipped
from janis_core.types import Array, Stdout
from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolArgument,
    ToolOutput,
    WildcardSelector,
    InputSelector,
    IndexOperator
)



# ------------- #
#  SECONDARIES  #
# ------------- #

# WORKFLOWS 

class SecondariesTestWF(Workflow):
    def id(self) -> str:
        return "SecondaryFileScatterTestWF"

    def friendly_name(self):
        return "WF which uses SecondaryFile types for workflow / tool IO"

    def constructor(self):
        self.input('inAlignments', BamBai())
        self.input('inAlignmentsArr', Array(BamBai()))
        self.input('inGunzipped', Gunzipped())

        # self.step(
        #     "stp1", 
        #     SecondariesTestTool(
        #         inp=self.inAlignments
        #     ), 
        # )
        # self.step(
        #     "stp2", 
        #     SecondariesTestTool(
        #         inp=self.stp1.out
        #     ), 
        # )
        # self.step(
        #     "stp3", 
        #     SecondariesReplacedTestTool(
        #         inp=self.inAlignments
        #     ), 
        # )
        self.step(
            "stp4", 
            SecondariesArrayTestTool(
                inp=self.inAlignmentsArr
            ), 
        )
        # self.step(
        #     "stp5", 
        #     GATKSplitReadsTestTool(
        #         bam=self.inAlignments
        #     ), 
        # )
        # self.step(
        #     "stp6", 
        #     NoSecondariesPresentAsTestTool(
        #         inp=self.inGunzipped
        #     ), 
        # )

        # self.output("outBamBai", source=self.stp1.out)
        # self.output("outBamBai2", source=self.stp2.out)



# TOOLS


class SecondariesTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp", 
                BamBai, 
                position=1),
        ]
    
    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument(
                InputSelector("inp"), 
                prefix='--inp', 
                position=2
            ),
            ToolArgument(
                IndexOperator(InputSelector("inp"), 0), 
                prefix='--inp-index-0', 
                position=3
            ),
            ToolArgument(
                IndexOperator(InputSelector("inp"), 1), 
                prefix='--inp-index-1', 
                position=4
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                BamBai,
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": ".bai"},
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class SecondariesReplacedTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesReplacedTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return []

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", BamBai, position=4)]

    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument("echo hello > out.bam", position=1, shell_quote=False),
            ToolArgument("&& echo there > out.bam.bai", position=2, shell_quote=False),
            ToolArgument("&& echo", position=3, shell_quote=False),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                BamBai,
                selector=WildcardSelector("*.bam"),
                secondaries_present_as={".bai": "^.bai"},
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

    

class SecondariesArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp", 
                Array(BamBai), 
                position=1),
        ]
    
    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument(
                InputSelector("inp"), 
                prefix='--inp', 
                position=2
            ),
            ToolArgument(
                IndexOperator(InputSelector("inp"), 0), 
                prefix='--inp-index-0', 
                position=3
            ),
            ToolArgument(
                IndexOperator(InputSelector("inp"), 1), 
                prefix='--inp-index-1', 
                position=4
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                'out',
                Stdout()
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class GATKSplitReadsTestTool(CommandTool):
    
    def friendly_name(self):
        return "TEST: GATKSplitReadsTestTool"

    def tool(self):
        return "GATKSplitReadsTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

    def inputs(self):
        return [
            ToolInput(
                "bam",
                BamBai,
                prefix="--input",
                position=1,
                secondaries_present_as={".bai": "^.bai"},
                doc="(-I:String) BAM/SAM/CRAM file containing reads  This argument must be specified at least once.",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                BamBai,
                glob=InputSelector("bam").basename(),
                doc="Bam",
                secondaries_present_as={".bai": "^.bai"},
            )
        ]


class NoSecondariesPresentAsTestTool(CommandTool):
    
    def friendly_name(self):
        return "TEST: NoSecondariesPresentAsTestTool"

    def tool(self):
        return "NoSecondariesPresentAsTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                Gunzipped(),
                position=8,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", FileTabix(), glob=InputSelector("inp")
            )
        ]