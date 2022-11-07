

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

from .types import SecondaryTestType


class SplitTextTestTool(PythonTool):
    @staticmethod
    def code_block(inp: str) -> dict:
        # list splits "abc" into ["a", "b", "c"]
        return {"out": list(inp)}

    def outputs(self):
        return [TOutput("out", Array(String()))]


class JoinArrayTestTool(PythonTool):
    @staticmethod
    def code_block(inp: Array(String())) -> dict:
        return {"out": " ".join(inp)}

    def outputs(self):
        return [TOutput("out", String())]


class SumTestTool(PythonTool):
    @staticmethod
    def code_block(inp1: int, inp2: int) -> dict:
        return {"out": inp1 + inp2}

    def outputs(self):
        return [TOutput("out", Int())]


class FileInputTestTool(PythonTool):
    @staticmethod
    def code_block(inp: File) -> dict:
        with open(inp) as f:
            content = f.read()

        return {"out": content}

    def outputs(self):
        return [TOutput("out", String())]


class SecondaryInputTestTool(PythonTool):
    @staticmethod
    def code_block(inp: SecondaryTestType) -> dict:
        with open(inp) as f:
            content = f.read()

        return {"out": content}

    def outputs(self):
        return [TOutput("out", String())]
