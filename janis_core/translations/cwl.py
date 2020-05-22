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
from io import StringIO

from typing import List, Dict, Optional, Tuple

import janis_core.utils.cwl_v1_0 as cwlgen
import ruamel.yaml

from janis_core.code.codetool import CodeTool
from janis_core.graph.steptaginput import full_lbl
from janis_core.tool.commandtool import CommandTool, ToolInput
from janis_core.tool.tool import Tool
from janis_core.translations.translationbase import TranslatorBase
from janis_core.types import (
    InputSelector,
    Selector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
    DataType,
    Directory,
)
from janis_core.types.common_data_types import Stdout, Stderr, Array, File, Filename
from janis_core.utils import first_value
from janis_core.utils.logger import Logger
from janis_core.utils.metadata import ToolMetadata

from janis_core.workflow.workflow import StepNode

CWL_VERSION = "v1.0"
SHEBANG = "#!/usr/bin/env cwl-runner"
yaml = ruamel.yaml.YAML()


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
            wf.id(),
            wf.friendly_name(),
            metadata.documentation,
            cwlVersion=CWL_VERSION,
            requirements=[],
        )

        w.inputs = [translate_input(i) for i in wf.input_nodes.values()]

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
                translate_step(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                )
            )

        w.outputs = [translate_output_node(o) for o in wf.output_nodes.values()]

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
                if i.value or i.default
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

        w.inputs = [translate_input(i) for i in wf.input_nodes.values()]

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
                translate_step(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                    use_run_ref=False,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
            )

        w.outputs = [translate_output_node(o) for o in wf.output_nodes]

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
        stdout = stdouts[0].stdoutname if len(stdouts) > 0 else None
        stderr = stderrs[0].stderrname if len(stderrs) > 0 else None

        if isinstance(stdout, InputSelector):
            stdout = translate_input_selector(stdout, code_environment=False)

        if isinstance(stderr, InputSelector):
            stderr = translate_input_selector(stderr, code_environment=False)

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
                    get_input_value_from_potential_selector_or_generator(
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
            translate_tool_input(i, inputsdict) for i in tool.inputs()
        )
        tool_cwl.outputs.extend(
            translate_tool_output(o, inputsdict=inputsdict, tool=tool.id())
            for o in tool.outputs()
        )

        args = tool.arguments()
        if args:
            tool_cwl.arguments.extend(
                translate_tool_argument(a) for a in tool.arguments()
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
                    coresMin="$(inputs.runtime_cpu ? inputs.runtime_cpu : 1)",
                    ramMin="$(inputs.runtime_memory ? Math.floor(1024 * inputs.runtime_memory) : 4096)",
                )
            )

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
            )
            for t in tool.inputs()
        )

        for output in tool.tool_outputs():
            if isinstance(output.outtype, Stdout):
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


def translate_input(inp):
    """

    :param inp:
    :type inp: janis_core.workflow.workflow.InputNode
    :return:
    """

    doc = inp.doc.doc if inp.doc else None

    return cwlgen.InputParameter(
        id=inp.id(),
        default=inp.default,
        secondaryFiles=inp.datatype.secondary_files(),
        format=None,
        streamable=None,
        doc=doc,
        inputBinding=None,
        label=None,
        type=inp.datatype.cwl_type(inp.default is not None),
    )


def translate_output_node(node):
    return translate_output(node, full_lbl(node.source[0], node.source[1]))


def translate_output(outp, source):
    ot = outp.datatype
    if isinstance(ot, Stdout):
        ot = ot.subtype or File()
    doc = outp.doc.doc if outp.doc else None

    return cwlgen.WorkflowOutputParameter(
        id=outp.id(),
        outputSource=source,
        secondaryFiles=outp.datatype.secondary_files(),
        streamable=None,
        doc=doc,
        type=ot.cwl_type(),
        linkMerge=None,
    )


def translate_tool_input(
    toolinput: ToolInput, inputsdict
) -> cwlgen.CommandInputParameter:

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
        value_from = get_input_value_from_potential_selector_or_generator(
            toolinput.default, code_environment=False, toolid=toolinput.id()
        )

    data_type = toolinput.input_type.cwl_type(default is not None)

    input_binding = cwlgen.CommandLineBinding(
        # load_contents=toolinput.load_contents,
        position=toolinput.position,
        prefix=toolinput.prefix,
        separate=toolinput.separate_value_from_prefix,
        itemSeparator=toolinput.separator,
        valueFrom=value_from,
        shellQuote=toolinput.shell_quote,
    )

    non_optional_dt_component = (
        [t for t in data_type if t != "null"][0]
        if isinstance(data_type, list)
        else data_type
    )

    # Binding array inputs onto the console
    # https://www.commonwl.org/user_guide/09-array-inputs/
    if isinstance(toolinput.input_type, Array) and isinstance(
        non_optional_dt_component, cwlgen.CommandInputArraySchema
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


def translate_tool_argument(argument):
    return cwlgen.CommandLineBinding(
        # load_contents=False,
        position=argument.position,
        prefix=argument.prefix,
        separate=argument.separate_value_from_prefix,
        # item_separator=None,
        valueFrom=get_input_value_from_potential_selector_or_generator(
            argument.value, code_environment=False
        ),
        shellQuote=argument.shell_quote,
    )


def translate_tool_output(output, inputsdict, **debugkwargs):

    doc = output.doc.doc if output.doc else None

    required_binding = not isinstance(output.output_type, (Stdout, Stderr))

    return cwlgen.CommandOutputParameter(
        id=output.tag,
        label=output.tag,
        secondaryFiles=prepare_tool_output_secondaries(output),
        # param_format=None,
        # streamable=None,
        doc=doc,
        outputBinding=cwlgen.CommandOutputBinding(
            glob=translate_to_cwl_glob(
                output.glob, inputsdict=inputsdict, outputtag=output.tag, **debugkwargs
            ),
            # load_contents=False,
            outputEval=prepare_tool_output_eval(output),
        )
        if required_binding
        else None,
        type=output.output_type.cwl_type(),
    )


def prepare_tool_output_eval(output):
    if not output.presents_as:
        return None
    return f"""${{
    self[0].basename="{output.presents_as}"
    return self
$}}
"""


def prepare_tool_output_secondaries(output):
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

    return f"""${{

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


def prepare_tool_input_secondaries(toolinp):
    if not toolinp.secondaries_present_as:
        return toolinp.input_type.secondary_files()

    secs = toolinp.secondaries_present_as
    tb = "    "
    formattedsecs = ",\n".join(
        f"""\
{4*tb}{{
{5*tb}location: resolveSecondary(self.location, "{secs.get(s, s)}"),
{5*tb}basename: resolveSecondary(self.basename, "{s}"),
{5*tb}class: "File",
{4*tb}}}"""
        for s in toolinp.input_type.secondary_files()
    )

    return f"""${{

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


def translate_step(
    step: StepNode,
    is_nested_tool=False,
    resource_overrides=Dict[str, str],
    use_run_ref=True,
    allow_empty_container=False,
    container_override=None,
):

    tool = step.tool
    if use_run_ref:
        prefix = "" if is_nested_tool else "tools/"
        run_ref = prefix + CwlTranslator.tool_filename(tool)
    else:
        from janis_core.workflow.workflow import Workflow

        has_resources_overrides = len(resource_overrides) > 0
        if isinstance(tool, Workflow):
            run_ref = CwlTranslator.translate_workflow_to_all_in_one(
                tool,
                with_resource_overrides=has_resources_overrides,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )
        elif isinstance(tool, CodeTool):
            run_ref = CwlTranslator.translate_code_tool_internal(
                tool,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )
        else:
            run_ref = CwlTranslator.translate_tool_internal(
                tool,
                True,
                with_resource_overrides=has_resources_overrides,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )

    cwlstep = cwlgen.WorkflowStep(
        id=step.id(),
        run=run_ref,
        # label=step.step.label,
        doc=step.doc.doc if step.doc else None,
        scatter=None,  # Filled by StepNode
        scatterMethod=None,  # Filled by StepNode
        in_=[],
        out=[],
    )

    cwlstep.out = [
        cwlgen.WorkflowStepOutput(id=o.tag) for o in step.tool.tool_outputs()
    ]

    ins = step.inputs()

    for k in ins:
        inp = ins[k]
        if k not in step.sources:
            if inp.intype.optional or inp.default:
                continue
            else:
                raise Exception(
                    f"Error when building connections for cwlstep '{step.id()}', "
                    f"could not find required connection: '{k}'"
                )

        edge = step.sources[k]
        ss = edge.slashed_source()
        link_merge = None

        if (
            ss is not None
            and not isinstance(ss, list)
            and isinstance(inp.intype, Array)
        ):
            start = edge.source().start
            outssval = start.outputs()
            source_type = (
                first_value(outssval)
                if len(outssval) == 1
                else outssval[edge.source().stag]
            ).outtype
            # has scattered = isinstance(start, StepNode) and start.scatter
            if not isinstance(source_type, Array) and not (
                isinstance(start, StepNode) and start.scatter
            ):
                ss = [ss]
                link_merge = "merge_nested"

        d = cwlgen.WorkflowStepInput(
            id=inp.tag,
            source=ss,
            linkMerge=link_merge,  # this will need to change when edges have multiple source_map
            valueFrom=None,
        )

        cwlstep.in_.append(d)

    for r in resource_overrides:
        cwlstep.in_.append(cwlgen.WorkflowStepInput(id=r, source=resource_overrides[r]))

    if step.scatter:
        if len(step.scatter.fields) > 1:
            cwlstep.scatterMethod = step.scatter.method.cwl()
        cwlstep.scatter = step.scatter.fields

    return cwlstep


## SELECTORS


def is_selector(selector):
    return issubclass(type(selector), Selector)


def get_input_value_from_potential_selector_or_generator(
    value, code_environment=True, **debugkwargs
):
    if value is None:
        return None
    if isinstance(value, str):
        return f'"{value}"' if code_environment else value
    elif isinstance(value, int) or isinstance(value, float):
        return value
    elif isinstance(value, Filename):
        # value.generated_filenamecwl() if code_environment else f"$({value.generated_filenamecwl()})"
        return (
            f'"{value.generated_filename()}"'
            if code_environment
            else value.generated_filename()
        )
    elif isinstance(value, StringFormatter):
        return translate_string_formatter(
            value, code_environment=code_environment, **debugkwargs
        )
    elif isinstance(value, InputSelector):
        return translate_input_selector(
            selector=value, code_environment=code_environment
        )
    elif isinstance(value, WildcardSelector):
        raise Exception(
            f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
        )
    elif isinstance(value, CpuSelector):
        return translate_cpu_selector(value)
    elif isinstance(value, MemorySelector):
        return translate_memory_selector(value)
    elif callable(getattr(value, "cwl", None)):
        return value.cwl()

    raise Exception("Could not detect type %s to convert to input value" % type(value))


def translate_input_selector(selector: InputSelector, code_environment):
    if not selector.input_to_select:
        raise Exception("No input was selected for input selector: " + str(selector))

    basename_extra = ".basename" if selector.use_basename else ""
    base = f"inputs.{selector.input_to_select}{basename_extra}"
    return base if code_environment else f"$({base})"


def translate_string_formatter(
    selector: StringFormatter, code_environment=True, **debugkwargs
):

    escapedFormat = selector._format.replace("\\", "\\\\")

    if len(selector.kwargs) == 0:
        return escapedFormat

    kwargreplacements = [
        f".replace(/{re.escape('{' +k + '}')}/g, {get_input_value_from_potential_selector_or_generator(v, code_environment=True, **debugkwargs)})"
        for k, v in selector.kwargs.items()
    ]
    return f'$("{escapedFormat}"' + "".join(kwargreplacements) + ")"


def translate_to_cwl_glob(glob, inputsdict, **debugkwargs):
    if not glob:
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
            else:
                expr = f"inputs.{glob.input_to_select}"
                if isinstance(intype, (File, Directory)):
                    expr = expr + ".basename"
                if tinp.default:
                    expr = f"(inputs.{glob.input_to_select} != null) ? {expr} : {tinp.default}"

                return f"$({expr})"

        return translate_input_selector(glob, code_environment=False)

    elif isinstance(glob, StringFormatter):
        return translate_string_formatter(glob)

    elif isinstance(glob, WildcardSelector):
        return glob.wildcard

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
    if not (inp and isinstance(inp, InputSelector)):
        return None

    if not inputsdict:
        raise Exception(
            f"Couldn't generate filename as an internal error occurred (inputsdict did not contain {inp.input_to_select})"
        )

    if inp.input_to_select not in inputsdict:
        raise Exception

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
