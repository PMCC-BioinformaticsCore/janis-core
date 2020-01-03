import argparse
from abc import ABC, abstractmethod
from textwrap import dedent
from typing import List

from janis_core.types import Boolean, Array, Int, Float, String

from janis_core.tool.tool import Tool, TOutput, TInput, ToolType, ToolTypes


class PythonTool(Tool, ABC):
    @abstractmethod
    def inputs(self) -> List[TInput]:
        pass

    @abstractmethod
    def outputs(self) -> List[TOutput]:
        pass

    def container(self):
        return "python:3.8.1"

    @staticmethod
    @abstractmethod
    def code_block(**kwargs):
        """
        This code block must be 100% self contained. All libraries and functions must be
        imported and declared from within this block.
        :param kwargs:
        :return:
        """
        pass

    # Other internal methods

    @classmethod
    def type(cls) -> ToolType:
        return ToolTypes.CodeTool

    def tool_inputs(self) -> List[TInput]:
        return self.inputs()

    def tool_outputs(self) -> List[TOutput]:
        return self.outputs()

    def generate_inputs_override(
        self, additional_inputs=None, with_resource_overrides=False, hints=None
    ):
        return {}

    def translate(
        self, translation: str, with_docker=True, with_resource_overrides=False
    ):
        pass

    def prepare_python_script(self):
        import inspect

        argkwargs = ", ".join(f"{t.id()}=args.{t.id()}" for t in self.inputs())

        codeblock_without_static = dedent(inspect.getsource(self.code_block)).split(
            "\n"
        )[1:]

        lines = [
            "import argparse, json",
            'cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")',
            *[self.generate_cli_binding_for_input(inp) for inp in self.tool_inputs()],
            *codeblock_without_static,
            "args = cli.parse_args()",
            "" f"result = code_block({argkwargs})",
            "",
            "print(json.dumps(result))",
        ]

        return "\n".join(lines)

    @staticmethod
    def generate_cli_binding_for_input(inp: TInput):
        params = [f'"--{inp.id()}"']
        intype = None

        required = not inp.intype.optional

        if isinstance(inp.intype, Int):
            intype = "int"
        elif isinstance(inp.intype, Float):
            intype = "float"
        elif isinstance(inp.intype, String):
            intype = "str"
        elif isinstance(inp.intype, Boolean):
            params.append("action='store_true'")
            required = False

        if intype:
            params.append("type=" + intype)

        if isinstance(inp.intype, Array):
            params.append("nargs='+'")

        if required:
            params.append("required=True")

        if inp.doc:
            escaped = inp.doc.replace("'", "\\'")
            params.append(f"help='{escaped}'")

        joined = ", ".join(params)
        return f"cli.add_argument({joined})"
