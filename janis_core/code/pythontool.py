from abc import ABC, abstractmethod
from textwrap import dedent
from typing import Dict, Any, Optional, Type

from janis_core.utils.docparser_info import parse_docstring
from janis_core.utils.metadata import ToolMetadata

from janis_core.code.codetool import CodeTool
from janis_core.tool.tool import TInput
from janis_core.types import (
    Boolean,
    Array,
    Int,
    Float,
    String,
    get_from_python_type,
    DataType,
    get_instantiated_type,
)

inspect_ignore_keys = {"self", "args", "kwargs", "cls", "template"}


class PythonTool(CodeTool, ABC):
    @staticmethod
    @abstractmethod
    def code_block(**kwargs) -> Dict[str, Any]:
        """
        This code block must be 100% self contained. All libraries and functions must be
        imported and declared from within this block.
        :param kwargs:
        :return:
        """
        pass

    def __init__(self, **connections):
        self._cached_input_signature = None
        super().__init__(metadata_class=ToolMetadata, **connections)

    # Other internal methods

    def inputs(self):
        if self._cached_input_signature is None:

            import inspect

            argspec = inspect.signature(self.code_block)
            docstrings = parse_docstring(self.code_block.__doc__)
            paramlist = docstrings.get("params", [])
            paramdocs = {
                p["name"]: p.get("doc").strip() for p in paramlist if "name" in p
            }

            unsupported_types = {}

            ins = []
            for inp in argspec.parameters.values():
                if inp.name in inspect_ignore_keys:
                    continue
                fdefault = inp.default
                optional = (fdefault is not inspect.Parameter.empty) or fdefault is None
                default = fdefault if optional else None

                defaulttype = type(fdefault) if fdefault is not None else None
                annotation = (
                    defaulttype
                    if inp.annotation is inspect.Parameter.empty
                    else inp.annotation
                )

                dt_type: Optional[DataType] = get_instantiated_type(
                    annotation, optional=optional
                )
                if not dt_type:
                    unsupported_types[inp.name] = annotation
                    continue

                ins.append(
                    TInput(
                        tag=inp.name,
                        intype=dt_type,
                        default=default,
                        doc=paramdocs.get(inp.name),
                    )
                )

            if unsupported_types:
                raise Exception(
                    f"Unsupported types for inputs: "
                    + ", ".join(f"{k}: {v}" for k, v in unsupported_types.items())
                )
            self._cached_input_signature = ins

        return self._cached_input_signature

    def base_command(self):
        return ["python", self.script_name()]

    def script_name(self):
        return f"{self.id()}-script.py"

    def container(self):
        return "python:3.8.1"

    def prepared_script(self):
        import inspect

        argkwargs = ", ".join(f"{t.id()}=args.{t.id()}" for t in self.inputs())

        codeblock_without_static = dedent(inspect.getsource(self.code_block)).split(
            "\n"
        )[1:]

        lines = [
            "import argparse, json",
            'cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")',
            *[self.generate_cli_binding_for_input(inp) for inp in self.tool_inputs()],
            "",
            *codeblock_without_static,
            "",
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
