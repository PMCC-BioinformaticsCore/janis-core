"""
CWL

This is one of the more complicated classes, it takes the janis in-memory representation of a tool,
and converts it into the equivalent CWL objects. Janis was built alongside testing for CWL, so a lot of
the concepts directly or pretty closely match. There are a few extra things that Janis has that need to
be mapped back.

This file is logically structured similar to the WDL equiv:

- Imports
- dump_cwl
- translate_workflow
- translate_tool (command tool)
- other translate methods
- selector helpers (InputSelector, WildcardSelector, CpuSelector, MemorySelector)
- helper methods
"""

## IMPORTS

import re
from typing import List, Dict, Optional, Any, Tuple, Union
from io import StringIO

from typing import List, Dict, Optional, Tuple

import cwl_utils.parser_v1_0 as cwlgen
import ruamel.yaml

from janis_core import ToolArgument
from janis_core.code.codetool import CodeTool
from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.operators.logical import IsDefined, If, RoundOperator
from janis_core.operators.standard import BasenameOperator, FirstOperator
from janis_core.tool.commandtool import CommandTool, ToolInput, ToolOutput
from janis_core.graph.steptaginput import full_lbl
from janis_core.tool.commandtool import CommandTool, ToolInput
from janis_core.tool.tool import Tool
from janis_core.translations.translationbase import TranslatorBase
from janis_core.operators import (
    InputSelector,
    Selector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
    Operator,
    InputNodeSelector,
    StepOutputSelector,
    TimeSelector,
    DiskSelector,
)
from janis_core.types.common_data_types import (
    Stdout,
    Stderr,
    Array,
    File,
    Filename,
    DataType,
    Directory,
)
from janis_core.utils import first_value
from janis_core.utils.logger import Logger
from janis_core.utils.metadata import ToolMetadata

from janis_core.workflow.workflow import StepNode, InputNode, OutputNode

CWL_VERSION = "v1.0"
SHEBANG = "#!/usr/bin/env cwl-runner"
yaml = ruamel.yaml.YAML()

STDOUT_NAME = "_stdout"
STDERR_NAME = "_stderr"


## TRANSLATION


class CwlTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="cwl")

        # class literal(str):
        #     pass
        #
        # def literal_presenter(dumper, data):
        #     return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        #
        # ruamel.yaml.add_representer(literal, literal_presenter)

    @staticmethod
    def walk_and_check_cwldict(c):
        return c.save()
        # primitives = (str, int, bool, float)
        # if isinstance(c, primitives):
        #     return c
        #
        # if isinstance(c, list):
        #     return [CwlTranslator.walk_and_check_cwldict(cc) for cc in c]
        #
        # elif isinstance(c, dict):
        #     return {k: CwlTranslator.walk_and_check_cwldict(v) for k, v in c.items()}
        # else:
        #     return CwlTranslator.walk_and_check_cwldict(c.save())

    @staticmethod
    def stringify_commentedmap(m):
        io = StringIO()
        yaml.dump(m, io)
        return io.getvalue()

    @staticmethod
    def stringify_translated_workflow(wf, format=True):
        saved = CwlTranslator.walk_and_check_cwldict(wf)
        formatted = SHEBANG + "\n" + CwlTranslator.stringify_commentedmap(saved)
        if format:
            from cwlformat.formatter import cwl_format

            formatted = cwl_format(formatted)

        return formatted

    @staticmethod
    def stringify_translated_tool(tool, format=True):
        saved = CwlTranslator.walk_and_check_cwldict(tool)

        formatted = SHEBANG + "\n" + CwlTranslator.stringify_commentedmap(saved)
        if format:
            from cwlformat.formatter import cwl_format

            formatted = cwl_format(formatted)

        return formatted

    @staticmethod
    def stringify_translated_inputs(inputs):
        return ruamel.yaml.dump(inputs, default_flow_style=False)

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return ["cwltool", "--validate", wfpath]

    @classmethod
    def translate_workflow(
        cls,
        wf,
        with_container=True,
        with_resource_overrides=False,
        is_nested_tool=False,
        is_packed=False,
        allow_empty_container=False,
        container_override=None,
    ) -> Tuple[cwlgen.Workflow, Dict[str, any]]:
        from janis_core.workflow.workflow import Workflow

        metadata = wf.metadata
        w = cwlgen.Workflow(
            id=wf.id(),
            label=wf.friendly_name(),
            doc=metadata.documentation,
            cwlVersion=CWL_VERSION,
            requirements=[],
            inputs=[],
            outputs=[],
            steps=[],
        )
        inputsdict = wf.inputs_map()

        w.inputs = [
            translate_workflow_input(i, inputsdict=inputsdict)
            for i in wf.input_nodes.values()
        ]

        resource_inputs = []
        if with_resource_overrides:
            resource_inputs = build_resource_override_maps_for_workflow(wf)
            w.inputs.extend(resource_inputs)

        w.steps = []

        for s in wf.step_nodes.values():
            resource_overrides = {}
            for r in resource_inputs:
                if not r.id.startswith(s.id()):
                    continue

                resource_overrides[r.id[(len(s.id()) + 1) :]] = r.id
            w.steps.append(
                translate_step_node(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                    allow_empty_container=allow_empty_container,
                )
            )

        w.outputs = [translate_workflow_output(o) for o in wf.output_nodes.values()]

        w.requirements.append(cwlgen.InlineJavascriptRequirement())
        w.requirements.append(cwlgen.StepInputExpressionRequirement())

        if wf.has_scatter:
            w.requirements.append(cwlgen.ScatterFeatureRequirement())
        if wf.has_subworkflow:
            w.requirements.append(cwlgen.SubworkflowFeatureRequirement())
        if wf.has_multiple_inputs:
            w.requirements.append(cwlgen.MultipleInputFeatureRequirement())

        tools = {}
        tools_to_build: Dict[str, Tool] = {
            s.tool.id(): s.tool for s in wf.step_nodes.values()
        }
        for t in tools_to_build:
            tool: Tool = tools_to_build[t]
            if isinstance(tool, Workflow):
                wf_cwl, subtools = cls.translate_workflow(
                    tool,
                    is_nested_tool=True,
                    with_container=with_container,
                    with_resource_overrides=with_resource_overrides,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
                tools[tool.versioned_id()] = wf_cwl
                tools.update(subtools)
            elif isinstance(tool, CommandTool):
                tool_cwl = cls.translate_tool_internal(
                    tool,
                    with_container=with_container,
                    with_resource_overrides=with_resource_overrides,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
                tools[tool.versioned_id()] = tool_cwl
            elif isinstance(tool, CodeTool):
                tool_cwl = cls.translate_code_tool_internal(
                    tool,
                    with_docker=with_container,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
                tools[tool.versioned_id()] = tool_cwl
            else:
                raise Exception(f"Unknown tool type: '{type(tool)}'")

        return w, tools

    @classmethod
    def build_inputs_file(
        cls,
        tool: Tool,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
    ) -> Dict[str, any]:
        from janis_core.workflow.workflow import Workflow

        ad = additional_inputs or {}
        values_provided_from_tool = {}
        if isinstance(tool, Workflow):
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value or (i.default and not isinstance(i.default, Selector))
            }

        inp = {
            i.id(): i.intype.cwl_input(
                ad.get(i.id(), values_provided_from_tool.get(i.id()))
            )
            for i in tool.tool_inputs()
            if i.default
            or not i.intype.optional
            or i.id() in ad
            or i.id() in values_provided_from_tool
        }

        if merge_resources:
            for k, v in cls.build_resources_input(
                tool, hints, max_cores, max_mem
            ).items():
                inp[k] = ad.get(k, v)

        return inp

    @classmethod
    def translate_workflow_to_all_in_one(
        cls,
        wf,
        with_resource_overrides=False,
        is_nested_tool=False,
        allow_empty_container=False,
        container_override=None,
    ) -> cwlgen.Workflow:

        metadata = wf.bind_metadata() or wf.metadata
        w = cwlgen.Workflow(
            wf.id(),
            wf.friendly_name(),
            metadata.documentation,
            cwlVersion=CWL_VERSION,
            requirements=[],
        )

        inputsdict = wf.inputs_map()

        w.inputs: List[cwlgen.InputParameter] = [
            translate_workflow_input(i, inputsdict=inputsdict)
            for i in wf.input_nodes.values()
        ]

        resource_inputs = []
        if with_resource_overrides:
            resource_inputs = build_resource_override_maps_for_workflow(wf)
            w.inputs.extend(resource_inputs)

        w.steps = []

        for s in wf.step_nodes.values():
            resource_overrides = {}
            for r in resource_inputs:
                if not r.id.startswith(s.id()):
                    continue

                resource_overrides[r.id[(len(s.id()) + 1) :]] = r.id

            w.steps.append(
                translate_step_node(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                    use_run_ref=False,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
            )

        w.outputs = [translate_workflow_output(o) for o in wf.output_nodes]

        w.requirements.append(cwlgen.InlineJavascriptRequirement())
        w.requirements.append(cwlgen.StepInputExpressionRequirement())

        if wf.has_scatter:
            w.requirements.append(cwlgen.ScatterFeatureRequirement())
        if wf.has_subworkflow:
            w.requirements.append(cwlgen.SubworkflowFeatureRequirement())
        if wf.has_multiple_inputs:
            w.requirements.append(cwlgen.MultipleInputFeatureRequirement())

        return w

    @classmethod
    def translate_tool_internal(
        cls,
        tool: CommandTool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override=None,
    ):
        metadata = tool.metadata if tool.metadata else ToolMetadata()

        stdout = STDOUT_NAME
        stderr = STDERR_NAME

        tool_cwl = cwlgen.CommandLineTool(
            id=tool.id(),
            baseCommand=tool.base_command(),
            label=tool.friendly_name() or tool.id(),
            doc=metadata.documentation,
            cwlVersion=CWL_VERSION,
            stdin=None,
            stderr=stderr,
            stdout=stdout,
            inputs=[],
            outputs=[],
            arguments=[],
            requirements=[],
        )

        # if any(not i.shell_quote for i in tool.inputs()):
        tool_cwl.requirements.append(cwlgen.ShellCommandRequirement())

        tool_cwl.requirements.extend([cwlgen.InlineJavascriptRequirement()])

        envs = tool.env_vars()
        if envs:
            lls = [
                cwlgen.EnvironmentDef(
                    k,
                    CwlTranslator.unwrap_expression(
                        value=v, code_environment=False, toolid=tool.id()
                    ),
                )
                for k, v in envs.items()
            ]
            tool_cwl.requirements.append(cwlgen.EnvVarRequirement(lls))

        inputs_that_require_localisation = [
            ti
            for ti in tool.inputs()
            if ti.localise_file
            and (
                isinstance(ti.input_type.received_type(), File)
                or (
                    issubclass(type(ti.input_type), Array)
                    and issubclass(type(ti.input_type.subtype()), File)
                )
            )
        ]
        if inputs_that_require_localisation:
            tool_cwl.requirements.append(
                cwlgen.InitialWorkDirRequirement(
                    [
                        cwlgen.Dirent(
                            entry="$(inputs.%s)" % ti.id(), entryname=ti.presents_as
                        )
                        for ti in inputs_that_require_localisation
                    ]
                )
            )

        if with_container:
            container = (
                CwlTranslator.get_container_override_for_tool(tool, container_override)
                or tool.container()
            )

            if container is not None:
                tool_cwl.requirements.append(
                    cwlgen.DockerRequirement(dockerPull=container)
                )
            elif not allow_empty_container:
                raise Exception(
                    f"The tool '{tool.id()}' did not have a container and no container override was specified. "
                    f"Although not recommended, Janis can export empty docker containers with the parameter "
                    f"'allow_empty_container=True' or --allow-empty-container"
                )

        inputsdict = {t.id(): t for t in tool.inputs()}
        tool_cwl.inputs.extend(
            translate_tool_input(i, inputsdict, tool) for i in tool.inputs()
        )
        tool_cwl.outputs.extend(
            translate_tool_output(o, inputsdict=inputsdict, tool=tool, toolid=tool.id())
            for o in tool.outputs()
        )

        args = tool.arguments()
        if args:
            tool_cwl.arguments.extend(
                translate_tool_argument(a, tool) for a in tool.arguments()
            )

        if with_resource_overrides:
            # work out whether (the tool of) s is a workflow or tool
            tool_cwl.inputs.extend(
                [
                    cwlgen.CommandInputParameter(id="runtime_memory", type="float?"),
                    cwlgen.CommandInputParameter(id="runtime_cpu", type="int?"),
                    # cwlgen.CommandInputParameter("runtime_disks", type="string?"),
                ]
            )

            tool_cwl.requirements.append(
                cwlgen.ResourceRequirement(
                    coresMin=CwlTranslator.unwrap_expression(
                        CpuSelector(), code_environment=False, tool=tool
                    ),
                    ramMin=CwlTranslator.unwrap_expression(
                        MemorySelector(), code_environment=False, tool=tool
                    ),
                    outdirMin=CwlTranslator.unwrap_expression(
                        DiskSelector(), code_environment=False, tool=tool
                    ),
                )
            )

            # Add ToolTimeLimit here

        return tool_cwl

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool: CodeTool,
        with_docker=True,
        allow_empty_container=False,
        container_override=None,
    ):

        stdouts = [
            o.outtype
            for o in tool.tool_outputs()
            if isinstance(o.outtype, Stdout) and o.outtype.stdoutname
        ]
        stderrs = [
            o.outtype
            for o in tool.tool_outputs()
            if isinstance(o.outtype, Stderr) and o.outtype.stderrname
        ]
        stdout = "python-capture.stdout"
        stderr = stderrs[0].stderrname if len(stderrs) > 0 else None

        scriptname = tool.script_name()
        inputsdict = {t.id(): t for t in tool.inputs()}

        if isinstance(stderr, InputSelector):
            stderr = translate_input_selector(stderr, code_environment=False)

        tool_cwl = cwlgen.CommandLineTool(
            id=tool.id(),
            baseCommand=tool.base_command(),
            label=tool.id(),
            doc="",  # metadata.documentation,
            cwlVersion=CWL_VERSION,
            stderr=stderr,
            stdout=stdout,
            inputs=[],
            outputs=[],
            requirements=[],
        )

        tool_cwl.inputs.extend(
            translate_tool_input(
                ToolInput(
                    t.id(),
                    input_type=t.intype,
                    prefix=f"--{t.id()}",
                    default=t.default,
                    doc=t.doc.doc if t.doc else None,
                ),
                inputsdict=inputsdict,
                tool=tool,
            )
            for t in tool.inputs()
        )

        for output in tool.tool_outputs():
            if isinstance(output.outtype, (Stdout, Stderr)):
                tool_cwl.outputs.append(
                    cwlgen.CommandOutputParameter(
                        id=output.tag, label=output.tag, type=output.outtype.cwl_type()
                    )
                )
                continue

            tool_cwl.outputs.append(
                cwlgen.CommandOutputParameter(
                    id=output.tag,
                    label=output.tag,
                    # param_format=None,
                    # streamable=None,
                    doc=output.doc.doc if output.doc else None,
                    outputBinding=cwlgen.CommandOutputBinding(
                        glob=stdout,
                        loadContents=True,
                        outputEval=cls.prepare_output_eval_for_python_codetool(
                            tag=output.tag, outtype=output.outtype
                        ),
                    ),
                    type=output.outtype.cwl_type(),
                )
            )

        tool_cwl.requirements.append(
            cwlgen.InitialWorkDirRequirement(
                listing=[
                    cwlgen.Dirent(entryname=scriptname, entry=tool.prepared_script())
                ]
            )
        )
        tool_cwl.requirements.append(cwlgen.InlineJavascriptRequirement())

        if with_docker:
            container = (
                CwlTranslator.get_container_override_for_tool(tool, container_override)
                or tool.container()
            )

            if container is not None:
                tool_cwl.requirements.append(
                    cwlgen.DockerRequirement(dockerPull=tool.container())
                )
            elif not allow_empty_container:
                raise Exception(
                    f"The tool '{tool.id()}' did not have a container. Although not recommended, "
                    f"Janis can export empty docker containers with the parameter 'allow_empty_container=True "
                    f"or --allow-empty-container"
                )

        return tool_cwl

    @staticmethod
    def prepare_output_eval_for_python_codetool(tag: str, outtype: DataType):

        requires_obj_capture = isinstance(outtype, (File, Directory))
        arraylayers = None
        if isinstance(outtype, Array) and isinstance(
            outtype.fundamental_type(), (File, Directory)
        ):
            requires_obj_capture = True
            base = outtype
            arraylayers = 0
            while isinstance(base, Array):
                arraylayers += 1
                base = outtype.subtype()

        out_capture = ""
        if requires_obj_capture:
            classtype = "File" if isinstance(outtype, File) else "Directory"
            fileout_generator = (
                lambda c: f"{{ class: '{classtype}', path: {c}, basename: {c}.substring({c}.lastIndexOf('/') + 1) }}"
            )

            if arraylayers:
                els = ["var els = [];"]

                base_var = f"v{arraylayers}"
                center = f"els.push({fileout_generator(base_var)};"

                def iteratively_wrap(center, iterable, layers_remaining):
                    var = f"v{layers_remaining}"
                    if layers_remaining > 1:
                        center = iteratively_wrap(center, var, layers_remaining - 1)
                    return f"for (var {var} of {iterable}) {{ {center} }}"

                out_capture = "\n".join(
                    [els, iteratively_wrap(center, "c", arraylayers)]
                )
            else:
                out_capture = fileout_generator("c")
        else:
            out_capture = "c"

        return f"""${{
var d = JSON.parse(self[0].contents)
if (!d) return null;
var c = d["{tag}"]
return {out_capture}
}}"""

    @classmethod
    def wrap_in_codeblock_if_required(cls, value, is_code_environment):
        return value if is_code_environment else f"$({value})"

    @classmethod
    def quote_values_if_code_environment(cls, value, is_code_environment):
        return f'"{value}"' if is_code_environment else value

    @classmethod
    def unwrap_selector_for_reference(cls, value):
        if value is None:
            return None
        if isinstance(value, str):
            return value
        elif isinstance(value, int) or isinstance(value, float):
            return value
        elif isinstance(value, InputNodeSelector):
            return value.id()
        elif isinstance(value, StepOutputSelector):
            return f"{value.node.id()}/{value.tag}"

        elif isinstance(value, InputSelector):
            return value.input_to_select

    @classmethod
    def unwrap_expression(
        cls,
        value,
        code_environment=True,
        selector_override=None,
        tool=None,
        for_output=False,
        **debugkwargs,
    ):
        if value is None:
            if code_environment:
                return "null"
            return None

        if isinstance(value, list):
            toolid = debugkwargs.get("tool_id", "unwrap_list_expression")
            inner = ", ".join(
                cls.unwrap_expression(
                    value[i],
                    code_environment=True,
                    selector_override=selector_override,
                    tool=tool,
                    tool_id=toolid + "." + str(i),
                )
                for i in range(len(value))
            )
            return cls.wrap_in_codeblock_if_required(
                f"[{inner}]", is_code_environment=code_environment
            )

        if isinstance(value, str):
            return CwlTranslator.quote_values_if_code_environment(
                value, code_environment
            )
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, Filename):
            # value.generated_filenamecwl() if code_environment else f"$({value.generated_filenamecwl()})"
            return CwlTranslator.quote_values_if_code_environment(
                value.generated_filename(), code_environment
            )

        elif isinstance(value, StringFormatter):
            return translate_string_formatter(
                value,
                selector_override=selector_override,
                code_environment=code_environment,
                tool=tool,
                **debugkwargs,
            )
        elif isinstance(value, InputNodeSelector):
            return translate_input_selector(
                InputSelector(value.id()),
                code_environment=code_environment,
                selector_override=selector_override,
            )
        elif isinstance(value, StepOutputSelector):
            raise Exception("Didn't expect to unwrap StepOutputSelector here")
            # sel = f"{value.node.id()}/{value.tag}"
            # if sel in selector_override:
            #     sel = selector_override[sel]
            # return CwlTranslator.wrap_in_codeblock_if_required(sel, code_environment)

        elif isinstance(value, MemorySelector):
            if not tool:
                raise Exception("Tool must be provided when unwrapping MemorySelector")
            toolmem = tool.memory({})

            if isinstance(toolmem, Operator) and any(
                isinstance(l, MemorySelector) for l in toolmem.get_leaves()
            ):
                raise Exception(
                    f"MemorySelector() should not be use used in tool.memory() for '{tool.id()}'"
                )
            ops = [InputSelector("runtime_memory")]
            if toolmem is not None:
                ops.append(toolmem)
            ops.append(4)
            # Note that 1GB = 953.674 MiB
            return cls.unwrap_expression(
                RoundOperator(953.674 * FirstOperator(ops)),
                code_environment=code_environment,
                tool=tool,
                **debugkwargs,
            )

        elif isinstance(value, CpuSelector):
            if not tool:
                raise Exception("Tool must be provided when unwrapping MemorySelector")
            toolcpu = tool.cpus({})

            if isinstance(toolcpu, Operator) and any(
                isinstance(l, CpuSelector) for l in toolcpu.get_leaves()
            ):
                raise Exception(
                    f"CpuSelector() should not be use used in tool.cpus() for '{tool.id()}'"
                )
            ops = [InputSelector("runtime_cpu")]
            if toolcpu is not None:
                ops.append(toolcpu)
            ops.append(1)
            return cls.unwrap_expression(
                FirstOperator(ops),
                code_environment=code_environment,
                tool=tool,
                **debugkwargs,
            )

        elif isinstance(value, TimeSelector):
            if not tool:
                raise Exception("Tool must be provided when unwrapping TimeSelector")
            tooltime = tool.time({})

            if isinstance(tooltime, Operator) and any(
                isinstance(l, TimeSelector) for l in tooltime.get_leaves()
            ):
                raise Exception(
                    f"TimeSelector() should not be use used in tool.time() for '{tool.id()}'"
                )
            ops = [InputSelector("runtime_seconds")]
            if tooltime is not None:
                ops.append(tooltime)
            ops.append(86400)
            return cls.unwrap_expression(
                FirstOperator(ops),
                code_environment=code_environment,
                tool=tool,
                **debugkwargs,
            )

        elif isinstance(value, DiskSelector):
            if not tool:
                raise Exception("Tool must be provided when unwrapping DiskSelector")
            tooltime = tool.disk({})

            if isinstance(tooltime, Operator) and any(
                isinstance(l, DiskSelector) for l in tooltime.get_leaves()
            ):
                raise Exception(
                    f"DiskSelector() should not be use used in tool.disk() for '{tool.id()}'"
                )
            ops = [InputSelector("runtime_disks")]
            if tooltime is not None:
                ops.append(tooltime)
            ops.append(20)
            return cls.unwrap_expression(
                FirstOperator(ops),
                code_environment=code_environment,
                tool=tool,
                **debugkwargs,
            )

        elif for_output and isinstance(value, (Stderr, Stdout)):
            # next few ones we rely on the globs being
            if isinstance(value, Stdout):
                return "self[0]"
            elif isinstance(value, Stderr):
                return "self[1]"

        elif isinstance(value, InputSelector):
            return translate_input_selector(
                selector=value,
                code_environment=code_environment,
                selector_override=selector_override,
            )
        elif isinstance(value, WildcardSelector):
            raise Exception(
                f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
            )
        elif isinstance(value, Operator):
            unwrap_expression_wrap = lambda exp: CwlTranslator.unwrap_expression(
                exp,
                code_environment=True,
                selector_override=selector_override,
                tool=tool,
                for_output=for_output,
                **debugkwargs,
            )
            return CwlTranslator.wrap_in_codeblock_if_required(
                value.to_cwl(unwrap_expression_wrap, *value.args),
                is_code_environment=code_environment,
            )
        elif callable(getattr(value, "cwl", None)):
            return value.cwl()
        # elif isinstance(value, Operator):

        raise Exception(
            "Could not detect type %s to convert to input value" % type(value)
        )

    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".cwl"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + "-inp.yml"

    @staticmethod
    def tool_filename(tool):
        return (tool.versioned_id() if isinstance(tool, Tool) else str(tool)) + ".cwl"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.yml"


def translate_workflow_input(inp: InputNode, inputsdict) -> cwlgen.InputParameter:
    """
    Translate a workflow InputNode into a cwlgen.InputParameter
    :param inp:
    :type inp: janis_core.workflow.workflow.InputNode
    :return:
    """

    doc = inp.doc.doc if inp.doc else None

    default = None
    dt = inp.datatype
    if inp.default is not None:
        if isinstance(inp.default, Selector):
            # CWL doesn't actually let us evaluate an expression here, even though it's
            # an inputBinding, we in fact will provide this as a default for the StepInput
            default = None
            dt = dt.copy()
            dt.optional = True
        else:
            default = inp.default

    return cwlgen.InputParameter(
        id=inp.id(),
        default=default,
        secondaryFiles=inp.datatype.secondary_files(),
        format=None,
        streamable=None,
        doc=doc,
        inputBinding=None,
        type=dt.cwl_type(default is not None),
    )


def translate_workflow_output(node: OutputNode) -> cwlgen.WorkflowOutputParameter:
    """
    Translate a workflow output node to a cwlgen.WorkflowOutputParameter
    :param node:
    :type node: OutputNode
    :return:
    """
    # we're going to need to transform this later to an operator
    source = CwlTranslator.unwrap_selector_for_reference(node.source)

    ot = node.datatype
    if isinstance(ot, Stdout):
        ot = ot.subtype or File()
    doc = node.doc.doc if node.doc else None

    return cwlgen.WorkflowOutputParameter(
        id=node.id(),
        outputSource=source,
        secondaryFiles=node.datatype.secondary_files(),
        streamable=None,
        doc=doc,
        type=ot.cwl_type(),
        linkMerge=None,
    )


def translate_tool_input(
    toolinput: ToolInput, inputsdict, tool
) -> cwlgen.CommandInputParameter:
    """
    Translate a ToolInput (commandtool / codetool) to a cwlgen.CommandInputParamaeter
    :param toolinput:
    :type toolinput: ToolInput
    :return:
    """

    default, value_from = toolinput.default, None
    intype = toolinput.input_type

    if isinstance(intype, Filename):

        if isinstance(intype.prefix, InputSelector):
            default = intype.generated_filename(
                inputs={intype.prefix.input_to_select: "generated"}
            )
            value_from = intype.generated_filename(
                inputs=prepare_filename_replacements_for(
                    intype.prefix, inputsdict=inputsdict
                )
            )
        else:
            default = intype.generated_filename()
    elif is_selector(default):
        default = None
        value_from = CwlTranslator.unwrap_expression(
            toolinput.default, code_environment=False, tool=tool, toolId=tool.id()
        )

    data_type = toolinput.input_type.cwl_type(default is not None)

    bind_to_commandline = toolinput.position is not None or toolinput.prefix is not None
    input_binding = (
        cwlgen.CommandLineBinding(
            # load_contents=toolinput.load_contents,
            position=toolinput.position,
            prefix=toolinput.prefix,
            separate=toolinput.separate_value_from_prefix,
            itemSeparator=toolinput.separator,
            valueFrom=value_from,
            shellQuote=toolinput.shell_quote,
        )
        if bind_to_commandline
        else None
    )

    non_optional_dt_component = (
        [t for t in data_type if t != "null"][0]
        if isinstance(data_type, list)
        else data_type
    )

    # Binding array inputs onto the console
    # https://www.commonwl.org/user_guide/09-array-inputs/
    if (
        bind_to_commandline
        and isinstance(toolinput.input_type, Array)
        and isinstance(non_optional_dt_component, cwlgen.CommandInputArraySchema)
    ):
        if toolinput.prefix_applies_to_all_elements:
            input_binding.prefix = None
            input_binding.separate = None
            nested_binding = cwlgen.CommandLineBinding(
                # load_contents=toolinput.load_contents,
                prefix=toolinput.prefix,
                separate=toolinput.separate_value_from_prefix,
                # item_separator=toolinput.item_separator,
                # value_from=toolinput.value_from,
                shellQuote=toolinput.shell_quote,
            )
            non_optional_dt_component.inputBinding = nested_binding

    doc = toolinput.doc.doc if toolinput.doc else None
    return cwlgen.CommandInputParameter(
        id=toolinput.tag,
        label=toolinput.tag,
        secondaryFiles=prepare_tool_input_secondaries(toolinput),
        # streamable=None,
        doc=doc,
        inputBinding=input_binding,
        default=default,
        type=data_type,
    )


def translate_tool_argument(argument: ToolArgument, tool) -> cwlgen.CommandLineBinding:
    """
    https://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineBinding

    :param argument: Tool argument to build command line for
    :return:
    """
    return cwlgen.CommandLineBinding(
        position=argument.position,
        prefix=argument.prefix,
        separate=argument.separate_value_from_prefix,
        valueFrom=CwlTranslator.unwrap_expression(
            argument.value, code_environment=False, tool=tool
        ),
        shellQuote=argument.shell_quote,
    )


def translate_tool_output(
    output: ToolOutput, inputsdict, tool, **debugkwargs
) -> cwlgen.CommandOutputParameter:
    """
    https://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputParameter

    :param output:
    :type output: ToolOutput
    :return:
    """

    doc = output.doc.doc if output.doc else None

    required_binding = not isinstance(output.output_type, (Stdout, Stderr))

    return cwlgen.CommandOutputParameter(
        id=output.tag,
        label=output.tag,
        secondaryFiles=prepare_tool_output_secondaries(output),
        doc=doc,
        outputBinding=prepare_tool_output_binding(
            output, inputsdict=inputsdict, tool=tool, **debugkwargs
        )
        if required_binding
        else None,
        type=output.output_type.cwl_type(),
    )


def requires_content(obj):
    if isinstance(obj, list):
        return any(requires_content(o) for o in obj)

    elif isinstance(obj, Operator) and obj.requires_contents():
        return True

    return False


def has_std(obj):
    if isinstance(obj, (Stderr, Stdout)):
        return True

    elif isinstance(obj, list):
        return any(has_std(o) for o in obj)

    elif isinstance(obj, Operator) and has_std(obj.get_leaves()):
        return True

    return False


def prepare_tool_output_binding(
    output: ToolOutput, inputsdict, tool, **debugkwargs
) -> cwlgen.CommandOutputBinding:

    loadcontents = requires_content(output.glob)
    requires_std = has_std(output.glob)

    return cwlgen.CommandOutputBinding(
        glob=[STDOUT_NAME, STDERR_NAME]
        if requires_std
        else translate_to_cwl_glob(
            output.glob, inputsdict, outputtag=output.tag, tool=tool, **debugkwargs
        ),
        outputEval=prepare_tool_output_eval(tool, output),
        loadContents=loadcontents,
    )


def prepare_tool_output_eval(tool, output: ToolOutput) -> Optional[str]:
    """
    Sometimes we want to rename it on the output, we do this by
    generating an output_eval expression and rewriting the basename

    :param output:
    :return:
    """

    if isinstance(output.glob, Operator):
        return CwlTranslator.unwrap_expression(
            output.glob, code_environment=False, tool=tool, for_output=True
        )

    if output.presents_as:
        return f"""\
${{
    self[0].basename="{output.presents_as}"
    return self
$}}
"""
    return None


def prepare_tool_output_secondaries(output) -> Optional[Union[str, List[str]]]:
    """
    Prepares the expressions / list of sec for a TOOL OUTPUT

    Secondary files CAN be pretty complicated.

    - If we're NOT 'secondaries_present_as', we return a
        simple list of secondary extension formats.

    - If we are using secondaries_present_as, we generate CWL expressions
        to write the basename with the correct format.

    :param output:
    :type output: ToolOutput
    :return:
    """

    if not output.secondaries_present_as:
        return output.output_type.secondary_files()

    secs = output.secondaries_present_as
    tb = "    "
    formattedsecs = ",\n".join(
        f"""\
{4*tb}{{
{5*tb}path: resolveSecondary(self.path, "{secs.get(s, s)}"),
{5*tb}basename: resolveSecondary(self.basename, "{s}"),
{5*tb}class: "File",
{4*tb}}}"""
        for s in output.output_type.secondary_files()
    )

    return [
        f"""${{

        function resolveSecondary(base, secPattern) {{
          if (secPattern[0] == "^") {{
            var spl = base.split(".");
            var endIndex = spl.length > 1 ? spl.length - 1 : 1;
            return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
          }}
          return base + secPattern
        }}
        return [
{formattedsecs}
        ];

}}"""
    ]


def prepare_tool_input_secondaries(inp: ToolInput) -> Optional[Union[str, List[str]]]:
    """
    Prepares the expressions / list of sec for a TOOL INPUT

    Secondary files CAN be pretty complicated.

    - If we're NOT 'secondaries_present_as', we return a
        simple list of secondary extension formats.

    - If we are using secondaries_present_as, we generate CWL expressions
        to write the basename with the correct format.

    :param inp:
    :type inp: ToolInput
    :return:
    """
    if not inp.secondaries_present_as:
        return inp.input_type.secondary_files()

    secs = inp.secondaries_present_as
    tb = "    "
    formattedsecs = ",\n".join(
        f"""\
{4*tb}{{
{5*tb}location: resolveSecondary(self.location, "{secs.get(s, s)}"),
{5*tb}basename: resolveSecondary(self.basename, "{s}"),
{5*tb}class: "File",
{4*tb}}}"""
        for s in inp.input_type.secondary_files()
    )

    return [
        f"""${{

        function resolveSecondary(base, secPattern) {{
          if (secPattern[0] == "^") {{
            var spl = base.split(".");
            var endIndex = spl.length > 1 ? spl.length - 1 : 1;
            return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
          }}
          return base + secPattern
        }}

        return [
{formattedsecs}
        ];

}}"""
    ]


def get_run_ref_from_subtool(
    tool: Tool,
    is_nested_tool,
    use_run_ref,
    resource_overrides=Dict[str, str],
    allow_empty_container=False,
    container_override=None,
):

    if use_run_ref:
        prefix = "" if is_nested_tool else "tools/"
        return prefix + CwlTranslator.tool_filename(tool)
    else:
        from janis_core.workflow.workflow import Workflow

        has_resources_overrides = len(resource_overrides) > 0
        if isinstance(tool, Workflow):
            return CwlTranslator.translate_workflow_to_all_in_one(
                tool,
                with_resource_overrides=has_resources_overrides,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )
        elif isinstance(tool, CodeTool):
            return CwlTranslator.translate_code_tool_internal(
                tool, allow_empty_container=allow_empty_container
            )
        else:
            return CwlTranslator.translate_tool_internal(
                tool,
                True,
                with_resource_overrides=has_resources_overrides,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )


def translate_step_node(
    step: StepNode,
    is_nested_tool=False,
    resource_overrides=Dict[str, str],
    use_run_ref=True,
    allow_empty_container=False,
    container_override=None,
) -> cwlgen.WorkflowStep:

    tool = step.tool

    # RUN REF
    run_ref = get_run_ref_from_subtool(
        tool,
        is_nested_tool=is_nested_tool,
        use_run_ref=use_run_ref,
        resource_overrides=resource_overrides,
        allow_empty_container=allow_empty_container,
        container_override=container_override,
    )

    # CONSTRUCTION

    cwlstep = cwlgen.WorkflowStep(
        id=step.id(),
        run=run_ref,
        label=tool.friendly_name(),
        doc=step.doc.doc if step.doc else None,
        scatter=None,  # Filled by StepNode
        scatterMethod=None,  # Filled by StepNode
        in_=[],
        out=[],
    )

    ## SCATTER

    if step.scatter:
        if len(step.scatter.fields) > 1:
            cwlstep.scatterMethod = step.scatter.method.cwl()
        cwlstep.scatter = step.scatter.fields

    ## OUTPUTS

    cwlstep.out = [
        cwlgen.WorkflowStepOutput(id=o.tag) for o in step.tool.tool_outputs()
    ]

    ## INPUTS

    for k, inp in step.inputs().items():
        if k not in step.sources:
            if inp.intype.optional or inp.default:
                continue
            else:
                raise Exception(
                    f"Error when building connections for cwlstep '{step.id()}', "
                    f"could not find required connection: '{k}'"
                )

        steptag_input: StepTagInput = step.sources[k]
        intype = inp.intype

        src: Edge = steptag_input.source()  # potentially single item or array
        ar_source = src if isinstance(src, list) else [src]
        has_multiple_sources = isinstance(src, list) and len(src) > 1
        array_input_from_single_source = False

        if len(ar_source) == 1:
            src = ar_source[0]

            ot = src.source.returntype()
            if (
                isinstance(intype, Array)
                and not isinstance(ot, Array)
                and not src.scatter
            ):
                array_input_from_single_source = True

        should_select_first_element = not (
            array_input_from_single_source or has_multiple_sources
        )

        has_operator = any(
            isinstance(src.source, Operator) for src in ar_source
        ) or any(
            (
                isinstance(src.source, InputNodeSelector)
                and isinstance(src.source.input_node.default, Selector)
            )
            for src in ar_source
        )

        source = None
        valuefrom = None
        link_merge = None
        default = None

        if not has_operator:
            unwrapped_sources: List[str] = []
            for stepinput in ar_source:
                src = stepinput.source
                unwrapped_sources.append(
                    CwlTranslator.unwrap_selector_for_reference(src)
                )

            source = (
                unwrapped_sources[0]
                if should_select_first_element
                else unwrapped_sources
            )
            if array_input_from_single_source:
                # Connect a solo value to an input that expects an array of that type
                # https://www.commonwl.org/user_guide/misc/
                link_merge = "merge_nested"

        else:
            prepare_alias = (
                lambda x: f"_{step.id()}_{k}_{re.sub('[^0-9a-zA-Z]+', '', x)}"
            )

            # two step process
            #   1. Look through and find ALL sources includng an operator's leaves
            #           Ensure these are connected using the alias
            #   2. Go through sources again, build up the expression

            param_aliasing = {}
            # Use a dict to ensure we don't double add inputs
            ins_to_connect: Dict[str, cwlgen.WorkflowStepInput] = {}

            processed_sources = []
            for i in range(len(ar_source)):
                stepinput = ar_source[i]
                src: any = stepinput.source
                if isinstance(src, InputNodeSelector) and isinstance(
                    src.input_node.default, Selector
                ):
                    src = If(IsDefined(src), src, src.input_node.default)

                processed_sources.append(src)

                if isinstance(src, Operator):
                    # we'll need to get the leaves and do extra mappings
                    for leaf in src.get_leaves():
                        if not isinstance(leaf, Selector):
                            # probably a python literal
                            continue
                        sel = CwlTranslator.unwrap_selector_for_reference(leaf)
                        alias = prepare_alias(sel)
                        param_aliasing[sel] = alias
                        ins_to_connect[alias] = cwlgen.WorkflowStepInput(
                            id=alias, source=sel
                        )
                else:
                    sel = CwlTranslator.unwrap_selector_for_reference(src)
                    alias = prepare_alias(sel)
                    param_aliasing[sel] = alias
                    ins_to_connect[alias] = cwlgen.WorkflowStepInput(
                        id=alias, source=sel
                    )

            cwlstep.in_.extend(ins_to_connect.values())

            # 2. Again

            src = (
                processed_sources[0]
                if should_select_first_element
                else [si for si in processed_sources]
            )
            valuefrom = CwlTranslator.unwrap_expression(
                src, code_environment=False, selector_override=param_aliasing, tool=tool
            )

        d = cwlgen.WorkflowStepInput(
            id=inp.tag,
            source=source,
            linkMerge=link_merge,  # this will need to change when edges have multiple source_map
            valueFrom=valuefrom,
            default=default,
        )

        cwlstep.in_.append(d)

    for r in resource_overrides:
        cwlstep.in_.append(cwlgen.WorkflowStepInput(id=r, source=resource_overrides[r]))

    return cwlstep


## SELECTORS


def is_selector(selector):
    return issubclass(type(selector), Selector)


def translate_input_selector(
    selector: InputSelector, code_environment, selector_override=None
):
    # TODO: Consider grabbing "path" of File

    sel = selector.input_to_select
    if not sel:
        raise Exception("No input was selected for input selector: " + str(selector))

    if selector_override and sel in selector_override:
        sel = selector_override[sel]

    basename_extra = ".basename" if selector.use_basename else ""
    base = f"inputs.{sel}{basename_extra}"
    return base if code_environment else f"$({base})"


def translate_string_formatter(
    selector: StringFormatter,
    selector_override,
    tool,
    code_environment=True,
    **debugkwargs,
):

    escapedFormat = selector._format.replace("\\", "\\\\")

    if len(selector.kwargs) == 0:
        return escapedFormat

    kwargreplacements = [
        f".replace(/{re.escape('{' +k + '}')}/g, {CwlTranslator.unwrap_expression(v, selector_override=selector_override, code_environment=True, tool=tool, **debugkwargs)})"
        for k, v in selector.kwargs.items()
    ]
    expr = f'"{escapedFormat}"' + "".join(kwargreplacements)
    if not code_environment:
        expr = f"$({expr})"
    return expr


def translate_to_cwl_glob(glob, inputsdict, tool, **debugkwargs):
    if glob is None:
        return None

    if not isinstance(glob, Selector):
        Logger.critical(
            "String globs are being phased out from tool output selections, please use the provided "
            "Selector (InputSelector or WildcardSelector) classes. " + str(debugkwargs)
        )
        return glob

    if isinstance(glob, InputSelector):

        if glob.input_to_select:
            if inputsdict is None or glob.input_to_select not in inputsdict:
                raise Exception(
                    "An internal error has occurred when generating the output glob for "
                    + glob.input_to_select
                )

            tinp: ToolInput = inputsdict[glob.input_to_select]
            intype = tinp.input_type
            if isinstance(intype, Filename):
                if isinstance(intype.prefix, InputSelector):
                    return intype.generated_filename(
                        inputs=prepare_filename_replacements_for(
                            intype.prefix, inputsdict=inputsdict
                        )
                    )
                else:
                    return intype.generated_filename()

            expr = glob
            if tinp.default:
                expr = If(IsDefined(glob), expr, tinp.default)

            return CwlTranslator.unwrap_expression(
                expr, code_environment=False, **debugkwargs
            )

    elif isinstance(glob, StringFormatter):
        return translate_string_formatter(glob, None, tool=tool)

    elif isinstance(glob, WildcardSelector):
        return glob.wildcard

    if isinstance(glob, Operator):
        # It's not terribly hard to do this, we'd have to change the output_eval
        # to use a combination of the presents_as AND
        if any(isinstance(o, WildcardSelector) for o in glob.get_leaves()):
            raise Exception(
                f"Janis does NOT currently support operations on the output glob for output '{glob}' in tool '{tool.id()}'"
            )
        return CwlTranslator.unwrap_expression(
            glob, code_environment=False, **debugkwargs
        )

    raise Exception("Unimplemented selector type: " + glob.__class__.__name__)


def translate_cpu_selector(selector: CpuSelector):
    return "$(inputs.runtime_cpu)"


def translate_memory_selector(selector: MemorySelector):
    return "$(Math.floor(inputs.runtime_memory))"


## OTHER HELPERS


def build_resource_override_maps_for_workflow(
    wf, prefix=None
) -> List[cwlgen.InputParameter]:
    from janis_core.workflow.workflow import Workflow

    # returns a list of key, value pairs
    inputs = []
    if not prefix:
        prefix = ""  # wf.id() + "."
    else:
        prefix += "_"

    for s in wf.step_nodes.values():
        tool: Tool = s.tool

        if isinstance(tool, CommandTool):
            tool_pre = prefix + s.id() + "_"
            inputs.extend(
                [
                    cwlgen.InputParameter(
                        id=tool_pre + "runtime_memory", type="float?"
                    ),
                    cwlgen.InputParameter(id=tool_pre + "runtime_cpu", type="int?"),
                    # cwlgen.InputParameter(tool_pre + "runtime_disks", type="string?"),
                ]
            )
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            inputs.extend(build_resource_override_maps_for_workflow(tool, tool_pre))

    return inputs


def prepare_filename_replacements_for(
    inp: Optional[InputSelector], inputsdict: Optional[Dict[str, ToolInput]]
) -> Optional[Dict[str, str]]:
    if inp is None or not isinstance(inp, InputSelector):
        return None

    if not inputsdict:
        raise Exception(
            f"Couldn't generate filename as an internal error occurred (inputsdict did not contain {inp.input_to_select})"
        )

    if inp.input_to_select not in inputsdict:
        raise Exception(
            f"The InputSelector '{inp.input_to_select}' did not select a valid input"
        )

    tinp = inputsdict.get(inp.input_to_select)
    intype = tinp.input_type

    if isinstance(intype, (File, Directory)):
        if isinstance(intype, File) and intype.extension:
            base = f'inputs.{tinp.id()}.basename.replace(/{intype.extension}$/, "")'
        else:
            base = f"inputs.{tinp.id()}.basename"
    else:
        base = "inputs." + tinp.id()

    if intype.optional:
        replacement = f'$(inputs.{tinp.id()} ? {base} : "generated")'
    else:
        replacement = f"$({base})"

    return {inp.input_to_select: replacement}
