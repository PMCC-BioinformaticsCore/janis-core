from itertools import groupby
from typing import Tuple, Dict, Union, List, Set

import inspect
from enum import Enum
from datetime import datetime, date

from janis_core.graph.steptaginput import StepTagInput

from janis_core.code.pythontool import PythonTool

from janis_core.types import DataType, UnionType

from janis_core.utils.logger import Logger

from janis_core.tool.commandtool import CommandTool, CommandToolBuilder
from janis_core.code.codetool import CodeTool
from janis_core.operators import Operator, InputSelector, StringFormatter, AliasSelector
from janis_core.workflow.workflow import (
    WorkflowBase,
    InputNode,
    OutputNode,
    InputNodeSelector,
    StepOutputSelector,
    StepNode,
)
from janis_core.tool.tool import Tool
from janis_core.translations import TranslatorBase


class HashOnlySet:
    """
    Keep adding anything hashable, but don't store duplicates if they're hash is the same
    """

    def __init__(self, d: dict = None):
        self.d = d or {}

    def add(self, el):
        return self.add_multiple([el])

    def add_multiple(self, els):
        added_c = 0
        for el in els:
            h = hash(el)
            if h in self.d:
                continue
            self.d[h] = el
            added_c += 1

        return added_c > 0

    def union(self, other_hashonly_set):
        return HashOnlySet({**self.d, **other_hashonly_set})

    def values(self):
        return list(self.d.values())

    def __contains__(self, item):
        return hash(item) in self.d

    def flatten_datatypes(self, dt: Union[DataType, List[DataType]]):
        """
        Pull out any subtypes from arrays, unions, etc
        """
        if not isinstance(dt, list):
            dt = [dt]

        flattened = []
        for d in dt:
            flattened.append(d)
            if d.is_array():
                flattened.extend(self.flatten_datatypes(d.subtype()))
            if isinstance(d, UnionType):
                flattened.extend(self.flatten_datatypes(d.subtypes))
        return flattened

    def consume_datatype(self, dt: Union[DataType, List[DataType]]):

        self.add_multiple(self.flatten_datatypes(dt))


class JanisTranslator(TranslatorBase):
    def __init__(self):
        super().__init__("janis-generator")

        self.import_tokens: HashOnlySet = HashOnlySet()

    @classmethod
    def get_class_name_from_tool(cls, tool: Tool):
        return tool.versioned_id().replace(".", "_").title()

    def translate_any_tool_internal(self, tool):
        if isinstance(tool, WorkflowBase):
            return self.generate_workflow_string(tool)
        elif isinstance(tool, CodeTool):
            return self.generate_code_tool_string(tool), {}
        elif isinstance(tool, CommandTool):
            return self.generate_command_tool_string(tool), {}
        raise NotImplementedError(
            f"Unrecognised tool type {tool}: {tool.__class__.__name__}"
        )

    def translate_workflow(
        self,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:
        str_wf, subtools = self.generate_workflow_string(workflow)

        nl = "\n"
        bigger_file = f"""
from datetime import datetime
from typing import List, Optional, Dict, Any

from janis_core import *
{self.prepare_imports()}

{nl.join(subtools.values())}

{str_wf}

{self.prepare_translate_command_at_bottom_of_generated_file(workflow)}
"""
        return bigger_file, {}

    def translate_tool_internal(
        self,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        tool_str = self.generate_command_tool_string(tool)

        bigger_file = f"""\
from datetime import datetime
from typing import List, Optional, Dict, Any

from janis_core import *
{self.prepare_imports()}

{tool_str}

{self.prepare_translate_command_at_bottom_of_generated_file(tool)}
"""

        return bigger_file

    def translate_code_tool_internal(
        self,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        tool_str = self.generate_code_tool_string(tool)

        bigger_file = f"""\
from datetime import datetime
from typing import List, Optional, Dict, Any

from janis_core import *
{self.prepare_imports()}

{tool_str}

{self.prepare_translate_command_at_bottom_of_generated_file(tool)}
"""
        return bigger_file

    def prepare_translate_command_at_bottom_of_generated_file(self, tool: Tool):
        extras = []
        if tool.has_tool_with_no_container():
            extras.append("allow_empty_container=True")

        extra_params = "".join(f", {param}" for param in extras)

        return f"""
if __name__ == "__main__":
    # or "cwl"
    {self.get_class_name_from_tool(tool)}().translate("wdl"{extra_params})
"""

    def prepare_imports(self):
        imports = []

        tokens_to_groupby = sorted(
            self.import_tokens.values(), key=lambda t: str(t.__module__)
        )

        for (imp, tokens) in groupby(tokens_to_groupby, lambda t: str(t.__module__)):

            imports.append(
                f"from {imp} import {', '.join(t.__class__.__name__ for t in tokens)}"
            )
        return "\n".join(imports)

    def generate_workflow_string(self, workflow) -> Tuple[str, Dict[str, str]]:
        w = self.get_class_name_from_tool(workflow)

        dts = [t.intype for t in workflow.tool_inputs()] + [
            t.outtype for t in workflow.tool_outputs()
        ]
        self.import_tokens.consume_datatype(dts)

        subtools = {}

        components = [self.convert_workflow_builder_initialiser(w, workflow)]

        for inp in workflow.input_nodes.values():
            components.append(self.convert_workflow_input(w, inp))

        for stp in workflow.step_nodes.values():
            if stp.tool.versioned_id() not in subtools:
                str_tool, more_subtools = self.translate_any_tool_internal(stp.tool)
                subtools = {
                    **more_subtools,
                    **subtools,
                    stp.tool.versioned_id(): str_tool,
                }

            components.append(self.convert_workflow_step(w, stp))

        for out in workflow.output_nodes.values():
            components.append(self.convert_workflow_output(w, out))

        return "\n".join(components), subtools

    def generate_command_tool_string(self, commandtool: CommandTool) -> str:
        if not isinstance(commandtool, CommandToolBuilder):
            commandtool = commandtool.to_command_tool_builder()

        # go through and add datatypes to to self.import_tokens
        dts = [t.input_type for t in commandtool.inputs()] + [
            t.output_type for t in commandtool.outputs()
        ]
        self.import_tokens.consume_datatype(dts)

        return f"{self.get_class_name_from_tool(commandtool)} = {self.convert_generic_class(commandtool)}"

    def generate_code_tool_string(self, tool: CodeTool):
        dts = [t.intype for t in tool.tool_inputs()] + [
            t.outtype for t in tool.tool_outputs()
        ]
        self.import_tokens.consume_datatype(dts)

        if isinstance(tool, PythonTool):
            return self.generate_python_tool_string(tool)
        else:
            NotImplementedError(
                f"Unsure how to generate janis string for code tool {tool}"
            )

    def generate_python_tool_string(self, tool: PythonTool):

        outs = self.get_string_repr(tool.outputs())

        cb = inspect.getsource(tool.code_block)

        return self.code_tool_format.format(
            tool_name=self.get_class_name_from_tool(tool),
            # input_fields=",".join(input_fields),
            code_block=cb,
            outputs=outs,
            identifier=tool.id(),
            version=tool.version(),
        )

    @classmethod
    def unwrap_expression(cls, expression):
        pass

    @staticmethod
    def stringify_translated_workflow(wf):
        try:
            import black

            try:
                return black.format_str(wf, mode=black.FileMode(line_length=82))
            except black.InvalidInput:
                Logger.warn(
                    "Check the generated Janis code carefully, as there might be a syntax error. You should report this error along with the workflow you're trying to generate from"
                )
        except ImportError:
            Logger.debug(
                "Janis can automatically format generated Janis code if you install black: https://github.com/psf/black"
            )

        return wf

        # return wf

    @staticmethod
    def stringify_translated_tool(tool):
        try:
            import black

            return black.format_str(tool, mode=black.FileMode(line_length=82))
        except ImportError:
            return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        return "inputs = " + str(inputs)

    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".py"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.versioned_id() + "-inputs.py"

    @staticmethod
    def tool_filename(tool):
        return tool.versioned_id() + ".py"

    @staticmethod
    def resources_filename(workflow):
        return workflow.versioned_id() + "-resources.py"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return ["python", wfpath]

    # Generic generator

    @staticmethod
    def get_string_repr(obj, workflow_id: str = None, get_string_repr_func=None):
        get_string_repr_func2 = lambda obj: (
            get_string_repr_func or JanisTranslator.get_string_repr
        )(obj, workflow_id=workflow_id, get_string_repr_func=get_string_repr_func)

        if isinstance(obj, list):
            inner = ", ".join(map(get_string_repr_func2, obj))
            return f"[{inner}]"
        elif isinstance(obj, dict):
            inner = ", ".join(
                f"{get_string_repr_func2(k)}: {get_string_repr_func2(v)}"
                for k, v in obj.items()
            )
            return f"{{{inner}}}"
        if isinstance(obj, str):
            nlreplaced = obj.replace("\n", "\\n").replace('"', "'")
            return f'"{nlreplaced}"'
        elif isinstance(obj, (int, float, bool, type(None))):
            return str(obj)
        elif isinstance(obj, datetime):
            return f"datetime({obj.year}, {obj.month}, {obj.day})"
        elif obj is None:
            return "None"
        elif isinstance(obj, Enum):
            return f"{obj.__class__.__name__}.{obj.name}"
        elif isinstance(obj, StringFormatter):
            inner_kwargs = "".join(
                f", {k}={get_string_repr_func2(v)}" for k, v in obj.kwargs.items()
            )
            return (
                f"StringFormatter({get_string_repr_func2(obj._format)}{inner_kwargs})"
            )
        elif isinstance(obj, Operator):
            return f"{obj.__class__.__name__}({', '.join(map(get_string_repr_func2, obj.args))})"
        elif isinstance(obj, AliasSelector):
            return f"{get_string_repr_func2(obj.inner_selector)}.as_type({get_string_repr_func2(obj.data_type)})"
        elif isinstance(obj, InputNodeSelector):
            return f"{workflow_id}.{obj.input_node.id()}"
        elif isinstance(obj, StepOutputSelector):
            return f"{workflow_id}.{obj.node.id()}.{obj.tag}"
        elif isinstance(obj, InputSelector) and workflow_id:
            return f"{workflow_id}.{obj.input_to_select}"
        elif isinstance(obj, date):
            return f"None"
        elif isinstance(obj, UnionType):
            return f"UnionType({', '.join(map(get_string_repr_func2, obj.subtypes))})"
        elif isinstance(
            obj, object
        ):  # any(isinstance(obj, T) for T in generic_convertible):
            return JanisTranslator.convert_generic_class(obj)

        return str(obj)

    @staticmethod
    def convert_generic_class(
        t, ignore_fields=None, get_string_repr_func=None, workflow_id: str = None
    ):
        options = []

        get_string_repr_func2 = lambda obj: (
            get_string_repr_func or JanisTranslator.get_string_repr
        )(obj, workflow_id=workflow_id)

        try:
            has_init_dict = not isinstance(
                t, (Tool, WorkflowBase, PythonTool, StepNode)
            ) and hasattr(t, "init_dictionary")
        except KeyError:
            has_init_dict = False

        if has_init_dict:
            options.extend(
                f"{k}={get_string_repr_func2(v)}"
                for k, v in t.init_dictionary().items()
            )
        else:
            ignore_fields = set(
                (ignore_fields if ignore_fields else []) + ["self", "args", "kwargs"]
            )

            params = inspect.signature(type(t).__init__).parameters
            param_map = {}
            if not isinstance(t, (StepNode, WorkflowBase)) and hasattr(
                t, "init_key_map"
            ):
                param_map = t.init_key_map
            # fields = fields_to_check if fields_to_check \
            #     else [f for f in dict(params).keys() if f not in ignore_fields]

            for fkey in params:
                if fkey in ignore_fields:
                    continue

                opts = params[fkey]

                t_key = param_map.get(fkey, fkey)
                if t_key is None:
                    continue

                if hasattr(t, t_key):
                    v = t.__getattribute__(t_key)
                else:
                    Logger.warn(
                        f"Object '{t.__class__.__name__}' didn't have attribute {t_key}, setting to None and it might get skipped"
                    )
                    v = None
                if (v is None and opts.default is None) or v == opts.default:
                    continue

                options.append(fkey + "=" + get_string_repr_func2(v))

        return f"{t.__class__.__name__}({', '.join(options)})"

    # Workflow translators

    def convert_workflow_builder_initialiser(
        self, workflow_identifier, workflow: WorkflowBase
    ):
        fields = [
            ("identifier", workflow.id()),
            ("friendly_name", workflow.friendly_name()),
            ("version", workflow.version()),
            ("tool_provider", workflow.tool_provider()),
            ("tool_module", workflow.tool_module()),
            ("doc", workflow.doc()),
        ]

        tb = 4 * " "
        mapped_fields = "\n".join(
            f"{tb}{k}={self.get_string_repr(v, workflow_id=workflow_identifier)},"
            for k, v in fields
            if v is not None
        )

        return f"""{workflow_identifier} = WorkflowBuilder(
{mapped_fields}
)
    """

    def convert_workflow_input(self, workflow_id, node: InputNode):

        fields = [
            ("default", node.default),
            ("value", node.value),
            ("doc", node.doc if node.doc.doc else None),
        ]
        tb = 4 * " "
        mapped_fields = "".join(
            f"\n{tb}{k}={self.get_string_repr(v, workflow_id=workflow_id)},"
            for k, v in fields
            if v is not None
        )

        return f"""{workflow_id}.input(
    "{node.id()}",
    {self.convert_generic_class(node.datatype)},{mapped_fields}
)
    """

    def convert_workflow_step(self, workflow_id, node: StepNode):
        extras = []
        if node.scatter is not None:
            extras.append(f"scatter={self.convert_generic_class(node.scatter)}")
        if node.when is not None:
            extras.append(
                f"when=${self.get_string_repr(node.when, workflow_id=workflow_id)}"
            )

        csp = 8 * " "
        connections = []
        for key, connection in node.sources.items():
            connection: StepTagInput = connection
            if (
                connection is None
                or isinstance(connection, (int, str, float, bool))
                or (
                    isinstance(connection, list)
                    and isinstance(connection[0], (int, str, float, bool))
                )
            ):
                connections.append(f"{csp}{key}={workflow_id}.{node.id()}_{key},")
                continue

            inp_source = (
                [e.source for e in connection.source_map]
                if (connection.multiple_inputs or len(connection.source_map) > 1)
                else connection.source_map[0].source
            )

            connections.append(
                f"{csp}{key}={self.get_string_repr(inp_source, workflow_id=workflow_id)},"
            )

        nl = "\n"
        return f"""
{workflow_id}.step(
    "{node.id()}",
    {self.get_class_name_from_tool(node.tool)}(
{nl.join(connections)}
    )
)
    """

    def convert_workflow_output(self, workflow_id, node: OutputNode):
        fields = [
            ("source", node.source),
            ("output_folder", node.output_folder),
            ("output_name", node.output_name),
            ("extension", node.extension),
            # ("doc", node.doc),
        ]
        tb = 4 * " "
        mapped_fields = "\n".join(
            f"{tb}{k}={self.get_string_repr(v, workflow_id=workflow_id)},"
            for k, v in fields
            if v is not None
        )
        return f"""{workflow_id}.output(
    "{node.id()}",
{mapped_fields}
)
    """

    # Code tool

    def prepare_type_annotation_for_code_tool(self, dt: DataType):
        inner = dt.name()
        if dt.is_array():
            inner = f"List[{self.prepare_type_annotation_for_code_tool(dt.subtype())}"
        if dt.optional:
            inner = f"Optional[{inner}]"
        return inner

    code_tool_format = """
class {tool_name}(PythonTool):
{code_block}

    def outputs(self) -> List[TOutput]:
        return {outputs}

    def id(self) -> str:
        return "{identifier}"

    def version(self):
        return "{version}"
    """


# extra templates


class ToolTemplateType(Enum):
    base = "regular"
    gatk4 = "gatk4"


def generate_gatk4_tooltemplatebase(gatk_command, inputs, outputs, metadata):
    io_prefix = " " * 12

    return gatk4_tool_template.format(
        gatk_command=gatk_command,
        inputs=",\n".join(io_prefix + s for s in inputs),
        outputs=",\n".join(io_prefix + s for s in outputs),
        metadata=metadata,
    )


def generate_regular_tooltemplatebase(
    toolname: str,
    name: str,
    friendly_name: str,
    tool_provider: str,
    base_command: Union[str, List[str]],
    inputs: List[str],  # =",\n".join((io_prefix + get_string_repr(i)) for i in ins),
    outputs: List[str],  # =",\n".join((io_prefix + get_string_repr(o)) for o in outs),
    metadata: str,  # =get_string_repr(metadata),
    version: str,
):
    import re

    io_prefix = " " * 12

    escapedversion = re.sub("[\W]", "_", str(version)) if version else str(None)

    return tool_template.format(
        toolname=toolname,
        name=name,
        friendly_name=friendly_name,
        tool_provider=tool_provider,
        base_command=base_command,
        inputs=",\n".join(io_prefix + s for s in inputs),
        outputs=",\n".join(io_prefix + s for s in outputs),
        metadata=metadata,
        version=version,
        escapedversion=escapedversion,
    )


tool_template = """
from abc import ABC
from datetime import datetime
from janis_core import (
    CommandTool, ToolInput, ToolOutput, File, Boolean, 
    String, Int, Double, Float, InputSelector, Filename, 
    ToolMetadata, InputDocumentation
)

class {name}Base(CommandTool, ABC):

    def friendly_name(self) -> str:
        return "{friendly_name}"

    def tool_provider(self):
        return "{tool_provider}"

    def tool(self) -> str:
        return "{toolname}"

    def base_command(self):
        return {base_command}

    def inputs(self):
        return [
{inputs}
        ]

    def outputs(self):
        return [
{outputs}
        ]

    def metadata(self):
        return {metadata}
        

class {name}_{escapedversion}({name}Base):
    def version(self):
        return "{version}"

    def container(self):
        return "{container}"
"""

gatk4_tool_template = """
from abc import ABC
from datetime import datetime
from janis_bioinformatics.tools.gatk4.gatk4toolbase import Gatk4ToolBase

from janis_core import (
    CommandTool, ToolInput, ToolOutput, File, Boolean, 
    String, Int, Double, Float, InputSelector, Filename, 
    ToolMetadata, InputDocumentation
)

class Gatk{gatk_command}Base(Gatk4ToolBase, ABC):

    @classmethod
    def gatk_command(cls):
        return "{gatk_command}"

    def friendly_name(self) -> str:
        return "GATK4: {gatk_command}"

    def tool(self) -> str:
        return "Gatk4{gatk_command}"

    def inputs(self):
        return [
{inputs}
        ]

    def outputs(self):
        return [
{outputs}
        ]

    def metadata(self):
        return {metadata}
"""
