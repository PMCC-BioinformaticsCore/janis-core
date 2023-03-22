

from typing import Optional

from janis_core import (
    ToolOutput,
    ToolInput,
    ToolArgument,
    CommandTool,
    Selector,
    Workflow,
    StringFormatter,
    InputSelector
)

from janis_core.types import (
    Stdout,
    File,
    String,
    Int
)

from ..testtools import MultiTypesInputPythonTool


class FilesDirectoriesToCreateTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            DirectoriesToCreateTestTool(inp=self.inStr)
        )
        self.step(
            "stp2", 
            FilesToCreateTestTool(
                inp1=self.inFile,
                inp2=self.inStr,
                inp3=self.inInt,
            )
        )
        self.step(
            "stp3", # code tool
            MultiTypesInputPythonTool(
                inp1=self.inFile,
                inp2=self.inStr,
                inp3=self.inInt,
            )
        )
        self.step(
            "stp4", 
            ExpressionToolTestTool(
                inp1=self.inFile,
                inp2=self.inStr,
                inp3=self.inInt,
            )
        )


    def friendly_name(self):
        return "TEST: FilesDirectoriesToCreateTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class DirectoriesToCreateTestTool(CommandTool):
    def tool(self) -> str:
        return "DirectoriesToCreateTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def directories_to_create(self) -> str | list[str]:
        return ['outdir']

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", str, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilesToCreateTestTool(CommandTool):
    def tool(self) -> str:
        return "FilesToCreateTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "bash"
    
    def files_to_create(self) -> dict[str, str | Selector]:
        return {
            'myscript.sh': 'this is a random script\n'
        }

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inp1", File, prefix='--inp1', position=1),
            ToolInput("inp2", String, prefix='--inp2', position=2),
            ToolInput("inp3", Int, prefix='--inp3', position=3),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    

class ExpressionToolTestTool(CommandTool):
    def tool(self) -> str:
        return "ExpressionToolTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "node check_value.js"
    
    def files_to_create(self) -> dict[str, str | Selector]:
        return {
            'check_value.js': '\nvar inputs = JSON.parse( process.argv[2] );\n\nvar ret = function(){\n  var value = 0;\n  if (inputs.number == 0 || inputs.number == 1) {\n    value = 2;\n  }\n  else {\n    value = inputs.number;\n  }\n  return {"out": value\n  }; }();\nprocess.stdout.write(JSON.stringify(ret));'
        }

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("script_file", File, presents_as='check_value.js'),
            ToolInput("inp1", File),
            ToolInput("inp2", String),
            ToolInput("inp3", Int),
        ]
    
    def arguments(self) -> Optional[list[ToolArgument]]:
        return [
            ToolArgument(
                StringFormatter(
                    format="'{{\"inp1\": {inp1}, \"inp2\": {inp2}, \"inp3\": {inp3}}}'", 
                    inp1=InputSelector('inp1'),
                    inp2=InputSelector('inp2'),
                    inp3=InputSelector('inp3'),
                )
            )
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


