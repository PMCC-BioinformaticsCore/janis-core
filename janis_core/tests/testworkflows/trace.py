

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
        self.input("inStringOpt", String(optional=True))
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
            StringFormatterTestTool(
                javaOptions=['--myname --ajeff'],
                compression_level=8
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


class StringFormatterTestTool(CommandTool):
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def tool(self) -> str:
        return "CatTestTool"

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
            # ToolInput("pg-tag", Boolean(optional=True), prefix="--add-output-sam-program-record",
            #           doc="If true, adds a PG tag to created SAM/BAM/CRAM files.")
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
        return [
            ToolOutput(
                "summaryMetrics",
                File(),
                glob=WildcardSelector(
                    "*.genotype_concordance_summary_metrics", select_first=True
                ),
            ),
            ToolOutput(
                "detailMetrics",
                File(),
                glob=WildcardSelector(
                    "*.genotype_concordance_detail_metrics", select_first=True
                ),
            ),
            ToolOutput(
                "contingencyMetrics",
                File(),
                glob=WildcardSelector(
                    "*.genotype_concordance_contingency_metrics", select_first=True
                ),
            ),
            # ToolOutput("vcf", VcfIdx(optional=True), glob=WildcardSelector("*.vcf"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


