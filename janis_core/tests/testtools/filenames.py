
from typing import Optional
from janis_core import CommandTool, ToolInput, ToolOutput, Stdout, File, Filename, InputSelector



# TODO outputs?????

class FilenameTestTool(CommandTool):
    def tool(self) -> str:
        return "FilenameTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

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

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilenameInputSelectorTestTool(CommandTool):
    def tool(self) -> str:
        return "FilenameInputSelectorTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

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

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"