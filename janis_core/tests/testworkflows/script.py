



from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector
)

from janis_core.types import (
    String,
    File,
    Array,
    Stdout,
    Filename,
    Boolean
)

from janis_bioinformatics.data_types.bam import BamBai


class ScriptTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inString', String())
        self.input('inSecondary', BamBai())
        
        self.input('inFileArr', Array(File()))
        self.input('inStringArr', Array(String()))
        self.input('inSecondaryArr', Array(BamBai()))

        self.step(
            "stp1", 
            ScriptTestTool1(
            )
        )

    def friendly_name(self):
        return "TEST: NamingTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class ScriptTestTool1(CommandTool):
    def tool(self) -> str:
        return "ScriptTestTool1"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inString", String()),        # note no position
            ToolInput("inFile", File, position=5),
            ToolInput(
                "inStringArrOpt",
                input_type=Array(String(), optional=True),
                prefix="--inStringArrOpt",
                prefix_applies_to_all_elements=True,
            ),
            ToolInput(
                "inStringFilename1",
                Filename(
                    prefix=InputSelector("inString"),
                    suffix="R1",
                    extension=".fastq.gz",
                ),
                prefix="--inStringFilename1",
            ),
            ToolInput(
                "inStringFilename2",
                Filename(
                    prefix=InputSelector("inString"),
                    suffix="R2",
                    extension=".fastq.gz",
                ),
                prefix="--inStringFilename2",
                doc="Write second read in a pair to FILE.",
            ),
            ToolInput(
                tag="inFlag",
                input_type=Boolean(optional=True),
                prefix="--inFlag",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

