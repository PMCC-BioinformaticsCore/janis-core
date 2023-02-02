

from janis_core import Workflow

from janis_core.types import (
    String,
    File,
)
from typing import Optional
from janis_core import CommandTool, ToolInput, ToolOutput, Stdout, File, Filename, InputSelector
from ..testtools import FilenameGeneratedTool



class FilenameTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inStrOpt', String(optional=True))
        self.input('inFileOpt', File(optional=True))

        self.step(
            "stp1", 
            FilenameTestTool(
                inp1=self.inFile,
                inp2=self.inStr,
            )
        )
        self.step(
            "stp2", 
            FilenameInputSelectorTestTool(
                inp1=self.inFile,
            )
        )
        self.step(
            "stp3", 
            FilenameGeneratedTool(
                inp=self.inStr,
                inpOptional=self.inStrOpt,
                fileInp=self.inFile,
                fileInpOptional=self.inFileOpt,
            )
        )
        self.step(
            "stp4", 
            FilenameCollectionTestTool(
                inp1=self.inFile,
            )
        )
        self.step(
            "stp5", 
            FilenameCollectionTestTool(
                inp1=self.inFile,
                inp4=self.inStr,
                inp5=self.inStr,
                inp6=self.inStr,
            )
        )

        self.output("outFilename", File, source=self.stp1.out)
        self.output("outFilenameInputSelector", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: FilenameTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    


class StdoutTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class FilenameTestTool(StdoutTestTool):
    def tool(self) -> str:
        return "FilenameTestTool"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp1", 
                File, 
                position=1
            ),
            ToolInput(
                "inp2", 
                Filename(), 
                position=2
            ),
        ]


class FilenameCollectionTestTool(CommandTool):
    
    def tool(self) -> str:
        return "FilenameCollectionTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp1", 
                File(), 
                position=1
            ),
            ToolInput(
                "inp3", 
                Filename(
                    prefix=InputSelector('inp1', remove_file_extension=True), 
                ), 
                position=2
            ),
            ToolInput(
                "inp4", 
                Filename(), 
                position=2
            ),
            ToolInput(
                "inp5", 
                Filename(extension='.csv'), 
                position=2
            ),
            ToolInput(
                "inp6", 
                Filename(suffix='.merged', extension='.csv'), 
                position=2
            ),
        ]

    def outputs(self):
        return [
            ToolOutput("out3", File(), glob=InputSelector("inp3") + ".csv"),    # ${inp1.simpleName + ".csv"}
            ToolOutput("out4", File(), glob=InputSelector("inp4") + ".csv"),   # ${(inp4 ? inp4 : "generated") + ".csv"}
            ToolOutput("out5", File(), glob=InputSelector("inp5")),   # ${(inp5 ? inp5 : "generated") + ".csv"}
            ToolOutput("out6", File(), glob=InputSelector("inp6")),   # ${inp6 + ".merged.csv"}
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class FilenameInputSelectorTestTool(StdoutTestTool):
    def tool(self) -> str:
        return "FilenameInputSelectorTestTool"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput(
                "inp1", 
                File, 
                position=1
            ),
            ToolInput(
                "inp2", 
                Filename(
                    prefix=InputSelector('inp1', remove_file_extension=True),
                    suffix='.processed',
                    extension='.txt'
                ), 
                position=2
            ),
        ]
