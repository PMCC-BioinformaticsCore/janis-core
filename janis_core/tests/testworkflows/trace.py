

from typing import Optional

from janis_core.types import File, Array, Stdout, String, Int, Filename
from janis_core import (
    Workflow, 
    ToolInput, 
    ToolArgument,
    ToolOutput, 
    CommandTool,
)

from janis_core import (
    If, 
    IsDefined,
    FirstOperator,
    WildcardSelector,
    InputSelector,
    StringFormatter,
    MemorySelector,
    JoinOperator
)


class EntityTraceTestWF(Workflow):

    def constructor(self):
        # self.input("inString", String(optional=True), value="")
        self.input("inFile", File())
        self.input("inFileArray", Array(File()))
        self.input("inString", String())
        self.input("inStringArrayOpt", Array(String(), optional=True))
        self.input("inStringOpt", String(optional=True))
        self.input("inIntOpt", Int(optional=True))
        self.input("inStringOptBackup", String(optional=True))
        self.input("inStringArray", Array(String()))

        # File -> File workflow input
        self.step(
            "stp1",
            CatTestTool(inp=self.inFile)
        )
        
        # File -> File step connection
        self.step(
            "stp2",
            CatTestTool(inp=self.stp1.out)
        )
        
        # If & IsDefined operators
        self.step(
            "stp3",
            EchoTestTool(inp=If(IsDefined(self.inStringOpt), self.inStringOpt, self.inStringOptBackup)),
        )
        
        # FirstOperator
        derived_string = FirstOperator(
            [
                self.inStringOpt,
                self.step(
                    "get_string",
                    EchoTestTool(inp="Some default value"),
                    when=self.inStringOpt.is_null(),
                ).out,
            ]
        )
        self.step(
            "stp4",
            EchoTestTool(inp=derived_string),
        )

        # IndexOperator
        self.step(
            "stp5", 
            EchoTestTool(inp=self.inStringArray[0])
        )

        self.step(
            "stp6",
            CatTestTool(inp=self.inFileArray),
            scatter='inp'
        )
        
        self.step(
            "stp7",
            CatTestTool(inp=self.stp6.out),
        )
        
        self.step(
            "stp8",
            InputUndefinedReferenceTestTool(
                # javaOptions=self.inStringArrayOpt,
                # compression_level=self.inIntOpt
            ),
        )
        self.step(
            "stp9",
            ArgumentUndefinedReferenceTestTool(
                # javaOptions=self.inStringArrayOpt,
                # compression_level=self.inIntOpt
            ),
        )
        self.step(
            "stp10",
            OutputUndefinedReferenceTestTool(
                # javaOptions=self.inStringArrayOpt,
                # compression_level=self.inIntOpt
            ),
        )

    def friendly_name(self):
        return "TEST: EntityTraceTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class EchoTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "EchoTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", String, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class CatTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "CatTestTool"

    def inputs(self) -> list[ToolInput]:
        return [ToolInput("inp", File, position=0)]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class InputUndefinedReferenceTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "InputUndefinedReferenceTestTool"

    def inputs(self):
        return [
            ToolInput("javaOptions", Array(String, optional=True)),
            ToolInput(
                "compression_level",
                Int(optional=True),
                doc="Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.",
            ),
            ToolInput(
                "javaOptionsRef",
                Filename(
                    prefix=InputSelector("javaOptions", remove_file_extension=True),
                    suffix=".fastq",
                    extension=".gz",
                ),
                prefix="--filenameRef",
                position=2,
            ),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"




class ArgumentUndefinedReferenceTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "ArgumentUndefinedReferenceTestTool"

    def inputs(self):
        return [
            ToolInput("javaOptions", Array(String, optional=True)),
            ToolInput(
                "compression_level",
                Int(optional=True),
                doc="Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.",
            ),
        ]

    def arguments(self):
        return [
            ToolArgument(
                StringFormatter(
                    "-Xmx{memory}G {compression} {otherargs}",
                    memory=MemorySelector() * 3 / 4,
                    compression=If(
                        IsDefined(InputSelector("compression_level")),
                        "-Dsamjdk.compress_level=" + InputSelector("compression_level"),
                        "",
                    ),
                    otherargs=JoinOperator(
                        FirstOperator([InputSelector("javaOptions"), []]), " "
                    ),
                ),
                prefix="--java-options",
                position=-1,
            )
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class OutputUndefinedReferenceTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "OutputUndefinedReferenceTestTool"

    def inputs(self):
        return [
            ToolInput("javaOptions", Array(String, optional=True)),
            ToolInput(
                "compression_level",
                Int(optional=True),
                doc="Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput("out", File(), glob=InputSelector("javaOptions")),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



