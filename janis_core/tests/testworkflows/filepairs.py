

from typing import Optional, Any
from janis_core.types import Array, File, Stdout
from janis_unix.data_types import ZipFile
from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector,
    IndexOperator
)




class FilePairsTestWF(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF"

    def friendly_name(self):
        return "TEST: FilePairsTestWF"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.input('inReadsArray', Array(FastqGzPairedEnd()))

        self.step(
            "stp1", 
            FilePairTestTool(
                reads=self.inReads
            ), 
        )
        self.step(
            "stp2", 
            FilePairTestTool(
                reads=self.inReadsArray
            ), 
            scatter='reads',
        )
        self.step(
            "stp3", 
            FilePairArrayTestTool(
                reads=self.inReadsArray
            ), 
        )


# TOOLS
class FilePairTestTool(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", FastqGzPairedEnd()),
            ToolInput(
                "read1",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 0),
                position=5,
            ),
            ToolInput(
                "read2",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 1),
                position=6,
            ),
        ]
    
    def outputs(self):
        return [
            ToolOutput(
                "out_R1",
                ZipFile(),
                selector=InputSelector("read1", remove_file_extension=True)
                + "_fastqc.zip",
            ),
            ToolOutput(
                "out_R2",
                ZipFile(),
                selector=InputSelector("read2", remove_file_extension=True)
                + "_fastqc.zip",
            ),

        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", Array(FastqGzPairedEnd())),
        ]
    
    def outputs(self):
        return [
            ToolOutput(
                "out",
                Stdout(),
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"




# redefining otherwise circular import error from janis_bioinformatics.
class FastqGz(File):
    def __init__(self, optional=False):
        super().__init__(
            optional=optional, extension=".fastq.gz", alternate_extensions={".fq.gz"}
        )

    @staticmethod
    def name():
        return "FastqGz"

    def doc(self):
        return (
            "FastqGz files are compressed sequence data with quality score, there are different types"
            "with no standard: https://en.wikipedia.org/wiki/FASTQ_format"
        )


# redefining otherwise circular import error from janis_bioinformatics.
class FastqGzPairedEnd(Array):
    def __init__(self, optional=False):
        super().__init__(FastqGz, optional=optional)

    @staticmethod
    def name():
        return "FastqGzPair"

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    def doc(self):
        return "Paired end FastqGz files"

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        super_is_valid = super().validate_value(meta, allow_null_if_not_optional)
        if not super_is_valid or meta is None:
            return super_is_valid

        return len(meta) == 2

    def invalid_value_hint(self, meta):
        prev = super().invalid_value_hint(meta)
        hints = []
        if prev:
            hints.append(prev)

        if meta is not None and len(meta) != 2:
            hints.append(f"There must be exactly 2 (found {len(meta)}) fastq files")
        return ", ".join(hints)

    @classmethod
    def ge(cls, file_paths: str, expected_sizes: list[int]):
        """

        :param file_paths: a string containing all file paths, separated by |
        :type file_paths: str
        :param expected_sizes: expected minimum sizes of all files
        :type expected_sizes: list[int]
        :return: a boolean value indicating if all files are bigger than or equal to their expected minimum sizes
        """
        files = file_paths.split("|")
        if len(files) != len(expected_sizes):
            return "Number of expected values don't match number of outputs"
        for file in files:
            unzipped = zipfile.ZipFile(file)
            if "R1" in unzipped.namelist()[0]:
                if os.path.getsize(file) < expected_sizes[0]:
                    return False
            else:
                if os.path.getsize(file) < expected_sizes[1]:
                    return False
        return True
