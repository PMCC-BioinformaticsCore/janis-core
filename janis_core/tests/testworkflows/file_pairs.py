

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


class FilePairsTestWF0(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF0"

    def friendly_name(self):
        return "TEST: FilePairsTestWF0"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairTestTool0(reads=self.inReads))


class FilePairsTestWF1(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF1"

    def friendly_name(self):
        return "TEST: FilePairsTestWF1"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairTestTool1(readsA=self.inReads, readsB=self.inReads))


class FilePairsTestWF2(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF2"

    def friendly_name(self):
        return "TEST: FilePairsTestWF2"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairTestTool2(reads=self.inReads))


class FilePairsTestWF3(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF2"

    def friendly_name(self):
        return "TEST: FilePairsTestWF2"

    def constructor(self):
        self.input('inFastq', FastqGz())
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairTestTool2(reads=self.inReads, read1=self.inFastq, read2=self.inFastq))
        

class FilePairsOptionalTestWF0(Workflow):
    def id(self) -> str:
        return "FilePairsOptionalTestWF0"

    def friendly_name(self):
        return "TEST: FilePairsOptionalTestWF0"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairOptionalTestTool0(reads=self.inReads))
        

class FilePairsOptionalTestWF1(Workflow):
    def id(self) -> str:
        return "FilePairsOptionalTestWF1"

    def friendly_name(self):
        return "TEST: FilePairsOptionalTestWF1"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairOptionalTestTool1(readsA=self.inReads, readsB=self.inReads))
        

class FilePairsOptionalTestWF2(Workflow):
    def id(self) -> str:
        return "FilePairsOptionalTestWF2"

    def friendly_name(self):
        return "TEST: FilePairsOptionalTestWF2"

    def constructor(self):
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairOptionalTestTool2(reads=self.inReads))
        

class FilePairsOptionalTestWF3(Workflow):
    def id(self) -> str:
        return "FilePairsOptionalTestWF3"

    def friendly_name(self):
        return "TEST: FilePairsOptionalTestWF3"

    def constructor(self):
        self.input('inFastq', FastqGz())
        self.input('inReads', FastqGzPairedEnd())
        self.step("stp1", FilePairOptionalTestTool2(reads=self.inReads, read1=self.inFastq, read2=self.inFastq))


class FilePairsArrayTestWF(Workflow):
    def id(self) -> str:
        return "FilePairsArrayTestWF"

    def friendly_name(self):
        return "TEST: FilePairsArrayTestWF"

    def constructor(self):
        self.input('inReadsArray', Array(FastqGzPairedEnd()))
        self.step("stp1", FilePairArrayTestTool1(read_pairs=self.inReadsArray))


class FilePairsArrayOptionalTestWF(Workflow):
    def id(self) -> str:
        return "FilePairsArrayOptionalTestWF"

    def friendly_name(self):
        return "TEST: FilePairsArrayOptionalTestWF"

    def constructor(self):
        self.input('inReadsArray', Array(FastqGzPairedEnd()))
        self.step("stp1", FilePairArrayOptionalTestTool1(read_pairs=self.inReadsArray))


class AllFilePairsTestWF(Workflow):
    def id(self) -> str:
        return "FilePairsTestWF"

    def friendly_name(self):
        return "TEST: FilePairsTestWF"

    def constructor(self):
        self.input('inFastq', FastqGz())
        self.input('inReads', FastqGzPairedEnd())
        self.input('inReadsOpt', FastqGzPairedEnd(optional=True))
        self.input('inReadsArray', Array(FastqGzPairedEnd()))
        self.input('inReadsArrayOpt', Array(FastqGzPairedEnd(), optional=True))

        self.step("stp1", FilePairTestTool1(reads=self.inReads))
        self.step("stp2", FilePairOptionalTestTool1(reads=self.inReadsOpt))
        self.step("stp3", FilePairArrayTestTool(read_pairs=self.inReadsArray))
        self.step("stp4", FilePairArrayOptionalTestTool(read_pairs=self.inReadsArrayOpt))


# TOOLS
class FilePairTestTool0(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool0"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", FastqGzPairedEnd(), position=1),
        ]
    
    def outputs(self):
        return [
            ToolOutput('stdout', Stdout())
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairTestTool1(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool1"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("readsA", FastqGzPairedEnd(), position=1, prefix='--prefix'),
            ToolInput("readsB", FastqGzPairedEnd(), position=1, prefix='--prefixeach', prefix_applies_to_all_elements=True),
        ]
    
    def outputs(self):
        return [
            ToolOutput('stdout', Stdout())
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairTestTool2(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool2"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", FastqGzPairedEnd(), position=1),
            ToolInput(
                "read1",
                FastqGz(optional=True),
                prefix='--reads-index-0',
                default=IndexOperator(InputSelector("reads"), 0),
                position=2,
            ),
            ToolInput(
                "read2",
                FastqGz(optional=True),
                prefix='--reads-index-1',
                default=IndexOperator(InputSelector("reads"), 1),
                position=3,
            ),
        ]
    
    def outputs(self):
        return [
            ToolOutput('stdout', Stdout()),
            ToolOutput(
                "out_R1",
                ZipFile(),
                selector=InputSelector("read1", remove_file_extension=True)
                + "_out.zip",
            ),
            ToolOutput(
                "out_R2",
                ZipFile(),
                selector=InputSelector("read2", remove_file_extension=True)
                + "_out.zip",
            )
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairOptionalTestTool0(CommandTool):
    def tool(self) -> str:
        return "FilePairOptionalTestTool0"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", FastqGzPairedEnd(optional=True), position=1),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    

class FilePairOptionalTestTool1(CommandTool):
    def tool(self) -> str:
        return "FilePairOptionalTestTool1"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("readsA", FastqGzPairedEnd(optional=True), position=1, prefix='--prefix'),
            ToolInput("readsB", FastqGzPairedEnd(optional=True), position=1, prefix='--prefixeach', prefix_applies_to_all_elements=True),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairOptionalTestTool2(CommandTool):
    def tool(self) -> str:
        return "FilePairOptionalTestTool2"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("reads", FastqGzPairedEnd(optional=True), position=1),
            ToolInput(
                "read1",
                FastqGz(optional=True),
                prefix='--reads-index-0',
                default=IndexOperator(InputSelector("reads"), 0),
                position=2,
            ),
            ToolInput(
                "read2",
                FastqGz(optional=True),
                prefix='--reads-index-1',
                default=IndexOperator(InputSelector("reads"), 1),
                position=3,
            ),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class FilePairArrayTestTool1(CommandTool):
    def tool(self) -> str:
        return "FilePairArrayTestTool1"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("read_pairs", Array(FastqGzPairedEnd()), position=1),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairArrayTestTool2(CommandTool):
    def tool(self) -> str:
        return "FilePairArrayTestTool2"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("read_pairs", Array(FastqGzPairedEnd()), position=1),
            ToolInput(
                "read_pairs0",
                FastqGz(optional=True),
                prefix='--read-pairs-0',
                default=IndexOperator(InputSelector("read_pairs"), 0),
                position=2,
            ),
            ToolInput(
                "read_pairs1",
                FastqGz(optional=True),
                prefix='--read-pairs-1',
                default=IndexOperator(InputSelector("read_pairs"), 1),
                position=3,
            ),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairArrayOptionalTestTool1(CommandTool):
    def tool(self) -> str:
        return "FilePairArrayOptionalTestTool1"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("read_pairs", Array(FastqGzPairedEnd(), optional=True), position=1),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class FilePairArrayOptionalTestTool2(CommandTool):
    def tool(self) -> str:
        return "FilePairArrayOptionalTestTool2"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("read_pairs", Array(FastqGzPairedEnd(), optional=True), position=1),
            ToolInput(
                "read_pairs0",
                FastqGz(optional=True),
                prefix='--read-pairs-0',
                default=IndexOperator(InputSelector("read_pairs"), 0),
                position=2,
            ),
            ToolInput(
                "read_pairs1",
                FastqGz(optional=True),
                prefix='--read-pairs-1',
                default=IndexOperator(InputSelector("read_pairs"), 1),
                position=3,
            ),
        ]
    
    def outputs(self):
        return [ToolOutput('stdout', Stdout())]

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
