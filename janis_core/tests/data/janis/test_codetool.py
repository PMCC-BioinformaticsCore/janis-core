

from janis_core import PythonTool
from janis_core import TOutput
from janis_core.types import File

class FileOutputPythonTestTool(PythonTool):
    @staticmethod
    def code_block(inp: File) -> dict:
        from shutil import copyfile

        copyfile(inp, "./out.file")

        return {"out": "./out.file"}

    def outputs(self):
        return [TOutput("out", File)]