


from typing import Optional
from janis_core.types import Array, File, Int
from janis_bioinformatics.data_types import FastqPairedEnd, BamBai

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector,
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
            ToolOutput("out_file", File(), selector=InputSelector("file")),
            ToolOutput("out_file_array", Array(File()), selector=InputSelector("file_array")),
            ToolOutput("out_file_optional", File(optional=True), selector=InputSelector("file_optional")),
            ToolOutput("out_file_array_optional", Array(File(), optional=True), selector=InputSelector("file_array_optional"))
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
            ToolOutput("out_filepair", FastqPairedEnd(), selector=InputSelector("filepair")),
            ToolOutput("out_filepair_array", Array(FastqPairedEnd()), selector=InputSelector("filepair_array")),
            ToolOutput("out_filepair_optional", FastqPairedEnd(optional=True), selector=InputSelector("filepair_optional")),
            ToolOutput("out_filepair_array_optional", Array(FastqPairedEnd(), optional=True), selector=InputSelector("filepair_array_optional"))
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
            ToolOutput("out_secondaries_array", Array(BamBai()), selector=InputSelector("secondaries_array")),
            ToolOutput("out_secondaries_optional", BamBai(optional=True), selector=InputSelector("secondaries_optional")),
            ToolOutput("out_secondaries_array_optional", Array(BamBai(), optional=True), selector=InputSelector("secondaries_array_optional"))
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"




# # redefining otherwise circular import error from janis_bioinformatics.
# class FastqGz(File):
#     def __init__(self, optional=False):
#         super().__init__(
#             optional=optional, extension=".fastq.gz", alternate_extensions={".fq.gz"}
#         )

#     @staticmethod
#     def name():
#         return "FastqGz"

#     def doc(self):
#         return (
#             "FastqGz files are compressed sequence data with quality score, there are different types"
#             "with no standard: https://en.wikipedia.org/wiki/FASTQ_format"
#         )


# # redefining otherwise circular import error from janis_bioinformatics.
# class FastqGzPairedEnd(Array):
#     def __init__(self, optional=False):
#         super().__init__(FastqGz, optional=optional)

#     @staticmethod
#     def name():
#         return "FastqGzPair"

#     def id(self):
#         if self.optional:
#             return f"Optional<{self.name()}>"
#         return self.name()

#     def doc(self):
#         return "Paired end FastqGz files"

#     def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
#         super_is_valid = super().validate_value(meta, allow_null_if_not_optional)
#         if not super_is_valid or meta is None:
#             return super_is_valid

#         return len(meta) == 2

#     def invalid_value_hint(self, meta):
#         prev = super().invalid_value_hint(meta)
#         hints = []
#         if prev:
#             hints.append(prev)

#         if meta is not None and len(meta) != 2:
#             hints.append(f"There must be exactly 2 (found {len(meta)}) fastq files")
#         return ", ".join(hints)

#     @classmethod
#     def ge(cls, file_paths: str, expected_sizes: list[int]):
#         """

#         :param file_paths: a string containing all file paths, separated by |
#         :type file_paths: str
#         :param expected_sizes: expected minimum sizes of all files
#         :type expected_sizes: list[int]
#         :return: a boolean value indicating if all files are bigger than or equal to their expected minimum sizes
#         """
#         files = file_paths.split("|")
#         if len(files) != len(expected_sizes):
#             return "Number of expected values don't match number of outputs"
#         for file in files:
#             unzipped = zipfile.ZipFile(file)
#             if "R1" in unzipped.namelist()[0]:
#                 if os.path.getsize(file) < expected_sizes[0]:
#                     return False
#             else:
#                 if os.path.getsize(file) < expected_sizes[1]:
#                     return False
#         return True
