

from typing import Optional
from janis_core.types import Array, Stdout

from janis_core.redefinitions.types import ZipFile, FastqGz, FastqGzPairedEnd
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


