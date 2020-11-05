from abc import ABC, abstractmethod
from textwrap import dedent, indent
from typing import Dict, Any, Optional, Type, List

from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.tool.documentation import InputDocumentation
from janis_core.types.data_types import NativeTypes

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
    DataType,
    File,
    Directory,
    get_instantiated_type,
    all_types,
)

inspect_ignore_keys = {"self", "args", "kwargs", "cls"}


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

            missing_annotations = set()
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
                if not annotation:
                    missing_annotations.add(inp.name)
                    continue

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
                        doc=InputDocumentation(paramdocs.get(inp.name)),
                    )
                )

            if missing_annotations:
                raise Exception(
                    f"The following types on the PythonTool '{self.id()}' were missing type annotations (REQUIRED): "
                    + ", ".join(missing_annotations)
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

    def prepared_script(self, translation: SupportedTranslation):
        import inspect

        nl = "\n"
        argkwargs = ", ".join(f"{t.id()}=args.{t.id()}" for t in self.inputs())

        codeblock_without_static = nl.join(
            dedent(inspect.getsource(self.code_block)).split(nl)[1:]
        )

        ins = self.tool_inputs()

        pt_decl = """\
class PythonTool:
    File = str
    Directory = str
"""

        type_annotation_declarations = nl.join(
            f"{k} = {v}"
            for k, v in {
                **{inp.intype.__class__.__name__: "str" for inp in ins},
                **{
                    t.__name__: NativeTypes.map_to_primitive(t.primitive()).__name__
                    for t in all_types
                },
                "Array": "List",
            }.items()
        )

        extra_param_parsing = ""

        if translation == SupportedTranslation.CWL:

            extra_param_parsing = """
from os import getcwd, path
cwd = getcwd()
def prepare_file_or_directory_type(file_or_directory, value):
    if value is None:
        return None
    if isinstance(value, list):
        return [prepare_file_or_directory_type(file_or_directory, v) for v in value]
    return {
        "class": file_or_directory,
        "path": path.join(cwd, value)
    }"""
            for out in self.outputs():
                st = (
                    out.outtype.fundamental_type()
                    if out.outtype.is_array()
                    else out.outtype
                )
                if not isinstance(st, (File, Directory)):
                    continue
                file_or_directory = "Directory" if isinstance(st, Directory) else "File"
                extra_param_parsing += f'\nresult["{out.id()}"] = prepare_file_or_directory_type("{file_or_directory}", result["{out.id()}"])'
            extra_param_parsing = indent(extra_param_parsing, 4 * " ")
        return f"""
import argparse, json, sys
from typing import Optional, List, Dict, Any
cli = argparse.ArgumentParser("Argument parser for Janis PythonTool")
{nl.join(self.generate_cli_binding_for_input(inp) for inp in ins)}

{type_annotation_declarations}
{pt_decl}


{codeblock_without_static}

try:
    args = cli.parse_args()
    result = code_block({argkwargs})
{extra_param_parsing}
    print(json.dumps(result))
except Exception as e:
    print(str(e), file=sys.stderr)
    raise
"""

    @staticmethod
    def generate_cli_binding_for_input(inp: TInput):
        params = [f'"--{inp.id()}"']

        intype = None

        bintype = inp.intype
        required = not inp.intype.optional

        if bintype.is_array():
            bintype = bintype.fundamental_type()
            params.append("nargs='+'")
            if required:
                params.append("default=[]")
                required = False

        if isinstance(bintype, Int):
            intype = "int"
        elif isinstance(bintype, Float):
            intype = "float"
        elif isinstance(bintype, (String, File, Directory)):
            intype = "str"
        elif isinstance(bintype, Boolean):
            params.append("action='store_true'")
            required = False

        if intype:
            params.append("type=" + intype)

        if required:
            params.append("required=True")

        if inp.doc and inp.doc.doc:
            escaped = inp.doc.doc.replace("'", "\\'").replace("\n", "\\n")
            params.append(f"help='{escaped}'")

        joined = ", ".join(params)
        return f"cli.add_argument({joined})"
