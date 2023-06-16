

from typing import Optional
from janis_core.redefinitions.types import BamBai, FileTabix, Gunzipped
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

        self.step(
            "stp1", 
            SecondariesTestTool(
                bam1=self.inAlignments,
                bam2=self.inAlignments,
                bam3=self.inAlignments,
            ), 
        )
        self.step(
            "stp2", 
            SecondariesOptionalTestTool(
                bam1=self.inAlignments,
            ), 
        )
        self.step(
            "stp3", 
            SecondariesArrayTestTool(
                bams1=self.inAlignmentsArr,
                bams2=self.inAlignmentsArr,
                bams3=self.inAlignmentsArr,
            ), 
        )
        self.step(
            "stp4", 
            SecondariesArrayOptionalTestTool(
                bams1=self.inAlignmentsArr,
            ), 
        )
        self.step(
            "stp5", 
            SecondariesReplacedTestTool(
                inp=self.inAlignments
            ), 
        )
        self.step(
            "stp6", 
            GATKSplitReadsTestTool(
                bam=self.inAlignments
            ), 
        )
        self.step(
            "stp7", 
            NoSecondariesPresentAsTestTool(
                inp=self.inGunzipped
            ), 
        )

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
                "bam1", 
                BamBai, 
                position=1
            ),
            ToolInput(
                "bam2", 
                BamBai, 
                position=2
            ),
            ToolInput(
                "bam3", 
                BamBai, 
                position=3
            ),
        ]
    
    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument(
                InputSelector("bam1"), 
                prefix='--arg1', 
                position=4
            ),
            ToolArgument(
                IndexOperator(InputSelector("bam1"), 0), 
                prefix='--arg2', 
                position=5
            ),
            ToolArgument(
                IndexOperator(InputSelector("bam1"), 1), 
                prefix='--arg3', 
                position=6
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


class SecondariesOptionalTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesOptionalTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "bam1", 
                BamBai(optional=True), 
                position=1
            )
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
                "bams1", 
                Array(BamBai), 
                position=1
            ),
            ToolInput(
                "bams2", 
                Array(BamBai), 
                prefix='--bams2',
                position=2
            ),
            ToolInput(
                "bams3", 
                Array(BamBai), 
                prefix='--bams3',
                prefix_applies_to_all_elements=True,
                position=3
            ),
        ]
    
    def arguments(self) -> list[ToolArgument]:
        return [
            ToolArgument(
                InputSelector("bams1"), 
                prefix='--bams-arg1', 
                position=4
            ),
            ToolArgument(
                IndexOperator(InputSelector("bams1"), 0), 
                prefix='--bams-arg2', 
                position=5
            ),
            ToolArgument(
                IndexOperator(InputSelector("bams1"), 1), 
                prefix='--bams-arg3', 
                position=6
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


class SecondariesArrayOptionalTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesArrayOptionalTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "bams1", 
                Array(BamBai, optional=True), 
                position=1
            )
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