

from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolArgument,
    ToolOutput,
    InputSelector,
    WildcardSelector,
    FirstOperator,
    IndexOperator
)
from janis_core.types import (
    Filename,
    File,
    String,
    Array,
)


from janis_unix import Csv, Tsv, TextFile

# WORKFLOW
class OutputCollectionTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inFileArray', Array(File))

        self.step("stp1", WildcardSelectorTestTool(inp=self.inFile))
        self.step("stp2", FilenameGenTestTool(inp=self.inFile))
        self.step("stp3", FilenameRefTestTool(inp=self.inFile))
        self.step("stp4", InputSelectorTestTool(inp=self.inFile, outputFilename='myfile.txt'))
        self.step("stp5", AddOperatorTestTool(inp=self.inFile))
        self.step("stp6", FilepairTestTool(inp=self.inFile))
        self.step("stp7", StringAddOperatorTestTool())
        self.step("stp8", MarkDuplicatesMetricsTestTool(outputPrefix='hello'))
        self.step("stp9", FilenameClashTestTool(reads=self.inFileArray ))

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


class FilepairTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: FilepairTestTool"

    def tool(self):
        return "FilepairTestTool"

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
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=[
                    (InputSelector("outputFilename") + '-R1.fastq'),
                    (InputSelector("outputFilename") + '-R2.fastq'),
                ]
            ),
        ]


class StringAddOperatorTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: StringAddOperatorTestTool"

    def tool(self):
        return "StringAddOperatorTestTool"

    def inputs(self):
        return [
            ToolInput(
                "outputPrefix",
                Filename(extension=".csv"),
                prefix="-o",
                doc="prefix of output summary csv",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", 
                Csv(),
                glob=InputSelector("outputPrefix") + ".csv"
            )
        ]


class MarkDuplicatesMetricsTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: MarkDuplicatesMetricsTestTool"

    def tool(self):
        return "MarkDuplicatesMetricsTestTool"

    def inputs(self):
        prefix = FirstOperator([InputSelector("outputPrefix"), "generated"])
        return [
            ToolInput("outputPrefix", String(optional=True)),
            ToolInput(
                "metricsFilename",
                Filename(prefix=prefix, suffix=".metrics", extension=".txt"),
                position=10,
                prefix="-M",
                doc="The output file to write marked records to.",
            ),
        ]

    def outputs(self):
        return [
            ToolOutput("metrics", Tsv(), glob=InputSelector("metricsFilename")),
        ]


class FilenameClashTestTool(CatToolBase):
    
    def friendly_name(self):
        return "TEST: FilenameClashTestTool"

    def tool(self):
        return "FilenameClashTestTool"

    def inputs(self):
        return [
            ToolInput("reads", Array(File)),
            ToolInput(
                "read1",
                File(optional=True),
                default=IndexOperator(InputSelector("reads"), 0),
                position=1,
            ),
            ToolInput(
                "read2",
                File(optional=True),
                default=IndexOperator(InputSelector("reads"), 1),
                position=2,
            ),
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out_R1_datafile",
                TextFile,
                selector=InputSelector("read1", remove_file_extension=True)
                + "_fastqc/fastqc_data.txt",
            ),
            ToolOutput(
                "out_R2_datafile",
                TextFile,
                selector=InputSelector("read2", remove_file_extension=True)
                + "_fastqc/fastqc_data.txt",
            ),
        ]


