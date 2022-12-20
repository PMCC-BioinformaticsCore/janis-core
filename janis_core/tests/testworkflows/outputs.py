

from typing import Optional

from janis_core import (
    CommandTool,
    ToolInput,
    ToolArgument,
    ToolOutput,
    InputSelector,
    WildcardSelector,
    Workflow
)
from janis_core.types import (
    Filename,
    File,
    String
)


# WORKFLOW
class OutputCollectionTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)

        self.step(
            "stp1", 
            WildcardSelectorTestTool(inp=self.inFile)
        )
        self.step(
            "stp2", 
            FilenameGenTestTool(inp=self.inFile)
        )
        self.step(
            "stp3", 
            FilenameRefTestTool(inp=self.inFile)
        )
        self.step(
            "stp4", 
            InputSelectorTestTool(
                inp=self.inFile,
                outputFilename='myfile.txt'
            )
        )
        
        self.step(
            "stp5", 
            AddOperatorTestTool(
                inp=self.inFile,
            )
        )

    def friendly_name(self):
        return "TEST: OutputCollectionTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# TOOLS

class CatToolBase(CommandTool):

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'


class WildcardSelectorTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: WildcardSelectorTestTool"

    def tool(self):
        return "WildcardSelectorTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                File(),
                position=1,
            ),
        ]
    
    def arguments(self):
        return [
            ToolArgument(
                "myfile.txt",
                prefix=">",
                position=2,
            )

        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=WildcardSelector("myfile.txt"),
            )
        ]


class FilenameGenTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: FilenameGenTestTool"

    def tool(self):
        return "FilenameGenTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                File(),
                position=1,
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    suffix=".recalibrated",
                    extension=".bam",
                ),
                prefix=">",
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=InputSelector("outputFilename"),
            )
        ]


class FilenameRefTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: FilenameRefTestTool"

    def tool(self):
        return "FilenameRefTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                File(),
                position=1,
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("inp", remove_file_extension=True),
                    suffix=".recalibrated",
                    extension=".bam",
                ),
                prefix=">",
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=InputSelector("outputFilename"),
            )
        ]


class InputSelectorTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: InputSelectorTestTool"

    def tool(self):
        return "InputSelectorTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                File(),
                position=1,
            ),
            ToolInput(
                "outputFilename",
                String(),
                prefix=">",
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=InputSelector("outputFilename"),
            )
        ]


class AddOperatorTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: AddOperatorTestTool"

    def tool(self):
        return "AddOperatorTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inp",
                File(),
                position=1,
            ),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("inp", remove_file_extension=True),
                ),
                prefix=">",
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=InputSelector("outputFilename")
                + ".gz",
            ),
        ]


