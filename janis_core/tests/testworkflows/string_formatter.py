

from typing import Optional, Any

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolArgument,
    InputSelector,
    FirstOperator,
    IndexOperator,
    MemorySelector,
    If,
    IsDefined,
    JoinOperator,
    StringFormatter
)

from janis_core.types import (
    Filename,
    File,
    String,
    Array,
    Int
)
from janis_core.redefinitions.types import BamBai


# WORKFLOW
class StringFormatterTestWF(Workflow):

    def constructor(self):
        self.input("inFile", File())
        self.input("inFileArr", Array(File()))
        self.input("inBamBai", BamBai())
        self.input("inBamBaiArr", Array(BamBai()))
        self.input("inStr", String())

        self.step(
            "stp1", 
            StringFormatterTestTool(
                javaOptions=['EXTRA_ARG1', 'EXTRA_ARG1'],
                compressionLevel=5,
            )
        )
        self.step(
            "stp2", 
            StringFormatterTestTool(
                compressionLevel=5,
            )
        )
        self.step(
            "stp3", 
            StringFormatterTestTool()
        )

    def friendly_name(self):
        return "TEST: UnwrapTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# TOOLS
class StringFormatterTestTool(CommandTool):

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

    def friendly_name(self):
        return "TEST: StringFormatterTestTool"

    def tool(self):
        return "StringFormatterTestTool"
    
    def cpus(self, hints: dict[str, Any]):
        return 4

    def memory(self, hints: dict[str, Any]):
        return 8

    def base_command(self) -> Optional[str | list[str]]:
        return 'echo'
    
    def inputs(self):
        return [
            ToolInput(
                "javaOptions", 
                Array(String, optional=True),
            ),
            ToolInput(
                "compressionLevel",
                Int(optional=True),
                prefix="--COMPRESSION_LEVEL",
                position=11,
                doc="Compression level for all compressed files created (e.g. BAM and GELI).",
            ),
        ]

    def arguments(self):
        return [
            ToolArgument(
                StringFormatter(
                    "-Xmx{memory}G {compression} {otherargs}",
                    memory=MemorySelector() * 3 / 4,
                    compression=If(
                        IsDefined(InputSelector("compressionLevel")),
                        "-Dsamjdk.compress_level=" + InputSelector("compressionLevel"),
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
        return []