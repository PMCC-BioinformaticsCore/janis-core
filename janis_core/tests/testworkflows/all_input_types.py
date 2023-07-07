


from typing import Optional
from janis_core.types import Array, File, Int
from janis_core.redefinitions.types import FastqPairedEnd, BamBai

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector,
    WildcardSelector
)

class AllInputTypesTestWF(Workflow):
    def id(self) -> str:
        return "AllInputTypesTestWF"

    def friendly_name(self):
        return "TEST: AllInputTypesTestWF"

    def constructor(self):
        self.input('in_file', File()) 
        self.input('in_file_array', Array(File()))
        self.input('in_file_optional', File(optional=True))
        self.input('in_file_array_optional', Array(File(),optional=True))
        self.input('in_nonfile', Int())
        self.input('in_nonfile_array', Array(Int()))
        self.input('in_nonfile_optional', Int(optional=True))
        self.input('in_nonfile_array_optional', Array(Int(), optional=True))
        self.input('in_nonfile_default', Int(), default=10)
        self.input('in_nonfile_array_default', Array(Int()), default=[1, 2, 3])
        self.input('in_filepair', FastqPairedEnd())
        self.input('in_filepair_array', Array(FastqPairedEnd()))
        self.input('in_filepair_optional', FastqPairedEnd(optional=True))
        self.input('in_filepair_array_optional', Array(FastqPairedEnd(), optional=True))
        self.input('in_secondaries', BamBai())
        self.input('in_secondaries_array', Array(BamBai()))
        self.input('in_secondaries_optional', BamBai(optional=True))
        self.input('in_secondaries_array_optional', Array(BamBai(), optional=True))

        self.step(
            "files", 
            FileTestTool(
                infile=self.in_file,
                infile_array=self.in_file_array,
                infile_optional=self.in_file_optional,
                infile_array_optional=self.in_file_array_optional,
            ), 
        )
        
        self.step(
            "nonfiles", 
            NonFileTestTool(
                nonfile=self.in_nonfile,
                nonfile_array=self.in_nonfile_array,
                nonfile_optional=self.in_nonfile_optional,
                nonfile_array_optional=self.in_nonfile_array_optional,
            ), 
        )
        
        self.step(
            "nonfiles_default", 
            NonFileDefaultTestTool(
                nonfile_default=self.in_nonfile_default,
                nonfile_array_default=self.in_nonfile_array_default,
            ), 
        )

        self.step(
            "filepairs", 
            FilePairTestTool(
                filepair=self.in_filepair,
                filepair_array=self.in_filepair_array,
                filepair_optional=self.in_filepair_optional,
                filepair_array_optional=self.in_filepair_array_optional,
            ), 
        )

        self.step(
            "secondaries", 
            SecondariesTestTool(
                secondaries=self.in_secondaries,
                secondaries_array=self.in_secondaries_array,
                secondaries_optional=self.in_secondaries_optional,
                secondaries_array_optional=self.in_secondaries_array_optional,
            ), 
        )




# TOOLS
class FileTestTool(CommandTool):
    def tool(self) -> str:
        return "FileTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('infile', File(), prefix='--file', position=1),
            ToolInput('infile_array', Array(File()), prefix='--file_array', position=2),
            ToolInput('infile_optional', File(optional=True), prefix='--file_optional', position=3),
            ToolInput('infile_array_optional', Array(File(), optional=True), prefix='--file_array_optional', position=4),
        ]
    
    def outputs(self):
        return [
            ToolOutput("out_file", File(), selector=InputSelector("infile")),
            ToolOutput("out_file_array", Array(File()), selector=InputSelector("infile_array")),
            ToolOutput("out_file_optional", File(optional=True), selector=InputSelector("infile_optional")),
            ToolOutput("out_file_array_optional", Array(File(), optional=True), selector=InputSelector("infile_array_optional"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class NonFileTestTool(CommandTool):
    def tool(self) -> str:
        return "NonFileTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('nonfile', Int(), prefix='--nonfile', position=1),
            ToolInput('nonfile_array', Array(Int()), prefix='--nonfile_array', position=2),
            ToolInput('nonfile_optional', Int(optional=True), prefix='--nonfile_optional', position=3),
            ToolInput('nonfile_array_optional', Array(Int(), optional=True), prefix='--nonfile_array_optional', position=4),
        ]
    
    def outputs(self):
        return [
            ToolOutput("out_nonfile", Int(), selector=InputSelector("nonfile")),
            ToolOutput("out_nonfile_array", Array(Int()), selector=InputSelector("nonfile_array")),
            ToolOutput("out_nonfile_optional", Int(optional=True), selector=InputSelector("nonfile_optional")),
            ToolOutput("out_nonfile_array_optional", Array(Int(), optional=True), selector=InputSelector("nonfile_array_optional"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class NonFileDefaultTestTool(CommandTool):
    def tool(self) -> str:
        return "NonFileDefaultTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('nonfile_default', Int(), prefix='--nonfile_default', position=1, default=1),
            ToolInput('nonfile_array_default', Array(Int()), prefix='--nonfile_array_default', position=2, default=[1, 2, 3]),
        ]
    
    def outputs(self):
        return [
            ToolOutput("out_nonfile_default", Int(), selector=InputSelector("nonfile_default")),
            ToolOutput("out_nonfile_array_default", Array(Int()), selector=InputSelector("nonfile_array_default")),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class FilePairTestTool(CommandTool):
    def tool(self) -> str:
        return "FilePairTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('filepair', FastqPairedEnd(), prefix='--filepair', position=1),
            ToolInput('filepair_array', Array(FastqPairedEnd()), prefix='--filepair_array', position=2),
            ToolInput('filepair_optional', FastqPairedEnd(optional=True), prefix='--filepair_optional', position=3),
            ToolInput('filepair_array_optional', Array(FastqPairedEnd(), optional=True), prefix='--filepair_array_optional', position=4),
        ]
    
    def outputs(self):
        return [
            ToolOutput("out_filepair", FastqPairedEnd(), selector=WildcardSelector("filepair")),
            ToolOutput("out_filepair_array", Array(FastqPairedEnd()), selector=WildcardSelector("filepair_array")),
            ToolOutput("out_filepair_optional", FastqPairedEnd(optional=True), selector=WildcardSelector("filepair_optional")),
            ToolOutput("out_filepair_array_optional", Array(FastqPairedEnd(), optional=True), selector=WildcardSelector("filepair_array_optional"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"



class SecondariesTestTool(CommandTool):
    def tool(self) -> str:
        return "SecondariesTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return ['echo']

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('secondaries', BamBai(), prefix='--secondaries', position=1),
            ToolInput('secondaries_array', Array(BamBai()), prefix='--secondaries_array', position=2),
            ToolInput('secondaries_optional', BamBai(optional=True), prefix='--secondaries_optional', position=3),
            ToolInput('secondaries_array_optional', Array(BamBai(), optional=True), prefix='--secondaries_array_optional', position=4),
        ]
    
    def outputs(self):
        return [
            ToolOutput("out_secondaries", BamBai(), selector=InputSelector("secondaries")),
            # ToolOutput("out_secondaries_array", Array(BamBai()), selector=InputSelector("secondaries_array")),
            ToolOutput("out_secondaries_optional", BamBai(optional=True), selector=InputSelector("secondaries_optional")),
            # ToolOutput("out_secondaries_array_optional", Array(BamBai(), optional=True), selector=InputSelector("secondaries_array_optional"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


