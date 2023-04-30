

from typing import Optional

from janis_core import (
    CommandTool,
    ToolInput,
    ToolArgument,
    InputSelector,
    FirstOperator,
    IndexOperator,
    Workflow
)
from janis_core.types import (
    Filename,
    File,
    String,
    Array
)

from janis_bioinformatics.data_types import BamBai


# TODO alias selector


# WORKFLOW
class UnwrapTestWF(Workflow):

    def constructor(self):
        self.input("inFile", File())
        self.input("inFileArr", Array(File()))
        self.input("inBamBai", BamBai())
        self.input("inBamBaiArr", Array(BamBai()))
        self.input("inStr", String())

        self.step(
            "stp1", 
            UnwrapTestTool(
                inFile=self.inFile,
                inFileArr=self.inFileArr,
                inBamBai=self.inBamBai,
                inBamBaiArr=self.inBamBaiArr,
                inStr=self.inStr,
            )
        )

    def friendly_name(self):
        return "TEST: UnwrapTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# TOOLS

class UnwrapTestTool(CommandTool):

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'

    def friendly_name(self):
        return "TEST: UnwrapTestTool"

    def tool(self):
        return "UnwrapTestTool"

    def inputs(self):
        return [
            ToolInput(
                "inFile",
                File(),
            ),
            ToolInput(
                "inFileArr",
                Array(File()),
            ),
            ToolInput(
                "inBamBai",
                BamBai(),
            ),
            ToolInput(
                "inBamBaiArr",
                Array(BamBai()),
            ),
            ToolInput(
                "inStr",
                String(),
            ),
            ToolInput(
                "filenameGen",
                Filename(extension=".gz"),
                prefix="--filenameGen",
                position=1,
            ),
            ToolInput(
                "filenameRef",
                Filename(
                    prefix=InputSelector("inFile", remove_file_extension=True),
                    suffix=".fastq",
                    extension=".gz",
                ),
                prefix="--filenameRef",
                position=2,
            ),
        ]
    
    def arguments(self):
        return [
            ToolArgument(
                InputSelector("inFile"),
                prefix="--inputSelectorProcess",
                position=3,
            ),
            ToolArgument(
                InputSelector("inStr"),
                prefix="--inputSelectorParam",
                position=4,
            ),
            ToolArgument(
                InputSelector("inFileArr"),
                prefix="--InputSelectorArray",
                position=8,
            ),
            ToolArgument(
                [1,2,3,4,5],
                prefix="--list",
                position=5,
            ),
            ToolArgument(
                InputSelector("inFile") + '.gz',
                prefix="--TwoValueOperator",
                position=6,
            ),
            ToolArgument(
                FirstOperator([InputSelector("inStr"), []]),
                prefix="--FirstOperator",
                position=7,
            ),
            ToolArgument(
                IndexOperator(InputSelector("inFileArr"), 0),
                prefix="--IndexOperatorArray",
                position=8,
            ),
            ToolArgument(
                IndexOperator(InputSelector("inBamBai"), 0),
                prefix="--IndexOperatorSecondariesBam",
                position=9,
            ),
            ToolArgument(
                IndexOperator(InputSelector("inBamBai"), 1),
                prefix="--IndexOperatorSecondariesBai",
                position=10,
            ),
            ToolArgument(
                IndexOperator(InputSelector("inBamBaiArr"), 0),
                prefix="--IndexOperatorArraySecondariesBams",
                position=11,
            ),
            ToolArgument(
                IndexOperator(InputSelector("inBamBaiArr"), 1),
                prefix="--IndexOperatorArraySecondariesBais",
                position=12,
            ),
        ]

    def outputs(self):
        return []
        #     ToolOutput(
        #         "out",
        #         File(),
        #         selector=WildcardSelector("myfile.txt"),
        #     )
        # ]
