

from janis_core import (
    PythonTool,
    TOutput
)

from janis_core.types import(
    Array, 
    String, 
    File,
    Int
) 

from .types import AppendedSecondaryTestType


class SplitTextPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: str) -> dict:
        # list splits "abc" into ["a", "b", "c"]
        return {"out": list(inp)}

    def outputs(self):
        return [TOutput("out", Array(String()))]


class JoinArrayPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: Array(String())) -> dict:
        return {"out": " ".join(inp)}

    def outputs(self):
        return [TOutput("out", String())]


class SumTestPythonTool(PythonTool):
    @staticmethod
    def code_block(inp1: int, inp2: int) -> dict:
        return {"out": inp1 + inp2}

    def outputs(self):
        return [TOutput("out", Int())]


class MultiTypesInputPythonTool(PythonTool):
    @staticmethod
    def code_block(inp1: File, inp2: String, inp3: Int) -> dict:
        from shutil import copyfile

        copyfile(inp1, "./out.file")

        return {"out": "./out.file"}

    def outputs(self):
        return [TOutput("out", File)]


# class ExpressionsPythonTool(PythonTool):
#     @staticmethod
#     def code_block(inp1: File, inp2: String, inp3: Int) -> dict:
#         from shutil import copyfile

#         copyfile(inp1, "./out.file")

#         return {"out": "./out.file"}
    
#     def arguments(self):
#         return [
#             ToolArgument(
#                 StringFormatter(
#                     "-Xmx{memory}G {compression} {otherargs}",
#                     memory=MemorySelector() * 3 / 4,
#                     compression=If(
#                         IsDefined(InputSelector("compression_level")),
#                         "-Dsamjdk.compress_level=" + InputSelector("compression_level"),
#                         "",
#                     ),
#                     otherargs=JoinOperator(
#                         FirstOperator([InputSelector("javaOptions"), []]), " "
#                     ),
#                 ),
#                 prefix="--java-options",
#                 position=-1,
#             )
#         ]

#     def outputs(self):
#         return [TOutput("out", File)]


class FileInputPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: File) -> dict:
        with open(inp) as f:
            content = f.read()

        return {"out": content}

    def outputs(self):
        return [TOutput("out", String())]


class FileOutputPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: File) -> dict:
        from shutil import copyfile

        copyfile(inp, "./out.file")

        return {"out": "./out.file"}

    def outputs(self):
        return [TOutput("out", File)]


class SecondaryInputPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: AppendedSecondaryTestType) -> dict:
        with open(inp) as f:
            content = f.read()

        return {"out": content}

    def outputs(self):
        return [TOutput("out", String())]
