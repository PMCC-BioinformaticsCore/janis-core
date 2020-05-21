"""
WDL

This is one of the more complicated classes, it takes the janis in-memory representation of a tool,
and converts it into the equivalent WDL objects. Python-wdlgen distances us from the memory representation
and the actual string-specification.

This file is logically structured similar to the cwl equiv:

- Imports
- dump_wdl
- translate_workflow
- translate_tool (command tool)
- other translate methods
- selector helpers (InputSelector, WildcardSelector, CpuSelector, MemorySelector)
- helper methods
"""

import json
from inspect import isclass
from typing import List, Dict, Any, Set, Tuple, Optional

import wdlgen as wdl

from janis_core.code.codetool import CodeTool
from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.tool.commandtool import CommandTool, ToolInput, ToolArgument, ToolOutput
from janis_core.tool.tool import Tool, TOutput
from janis_core.translations.translationbase import TranslatorBase
from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    StringFormatter,
    String,
    Selector,
    Directory,
)
from janis_core.types.common_data_types import (
    Stdout,
    Stderr,
    Array,
    Boolean,
    Filename,
    File,
)
from janis_core.utils import first_value, recursive_2param_wrap, find_duplicates
from janis_core.utils.generators import generate_new_id_from
from janis_core.utils.logger import Logger
from janis_core.utils.scatter import ScatterDescription, ScatterMethods
from janis_core.utils.validators import Validators
from janis_core.utils.secondary import (
    split_secondary_file_carats,
    apply_secondary_file_format_to_filename,
)

# from janis_core.tool.step import StepNode


## PRIMARY TRANSLATION METHODS

SED_REMOVE_EXTENSION = "| sed 's/\\.[^.]*$//'"
REMOVE_EXTENSION = (
    lambda x, iterations: f"`echo '{x}' {iterations * SED_REMOVE_EXTENSION}`"
)


class CustomGlob(Selector):
    def __init__(self, expression):
        self.expression = expression


class WdlTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="wdl")

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf.get_string()

    @staticmethod
    def stringify_translated_tool(tool):
        return tool.get_string()

    @staticmethod
    def stringify_translated_inputs(inputs):
        return json.dumps(inputs, sort_keys=True, indent=4, separators=(",", ": "))

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return ["java", "-jar", "$womtooljar", "validate", wfpath]

    @classmethod
    def translate_workflow(
        cls,
        wfi,
        with_container=True,
        with_resource_overrides=False,
        is_nested_tool=False,
        allow_empty_container=False,
        container_override=None,
    ) -> Tuple[any, Dict[str, any]]:
        """
        Translate the workflow into wdlgen classes!


        :param with_resource_overrides:
        :param with_container:
        :param is_nested_tool:
        :return:
        """
        # Import needs to be here, otherwise we end up circularly importing everything
        # I need the workflow for type comparison
        from janis_core.workflow.workflow import Workflow

        wf: Workflow = wfi

        # Notes:
        #       All wdlgen classes have a .get_string(**kwargs) function
        #       The wdlgen Workflow class requires a

        # As of 2019-04-16: we use development (Biscayne) features
        # like Directories and wdlgen uses the new input {} syntax
        w = wdl.Workflow(wf.id(), version="development")
        tools: List[Tool] = [s.tool for s in wf.step_nodes.values()]

        inputs = list(wf.input_nodes.values())
        steps = list(wf.step_nodes.values())
        outputs = list(wf.output_nodes.values())

        wtools = {}  # Store all the tools by their name in this dictionary
        tool_aliases, step_aliases = build_aliases(
            wf.step_nodes.values()
        )  # Generate call and import aliases

        # Convert self._inputs -> wdl.Input
        for i in inputs:
            dt = i.datatype
            expr = None
            if isinstance(i.datatype, Filename):
                # expr = f'"{i.datatype.generated_filename()}"'
                dt = String(optional=True)

            wd = dt.wdl(has_default=i.default is not None)

            w.inputs.append(
                wdl.Input(
                    data_type=wd, name=i.id(), expression=expr, requires_quotes=False
                )
            )

            is_array = isinstance(i.datatype, Array)
            if i.datatype.secondary_files() or (
                is_array and i.datatype.subtype().secondary_files()
            ):
                secs = (
                    i.datatype.secondary_files()
                    if not is_array
                    else i.datatype.subtype().secondary_files()
                )

                w.inputs.extend(
                    wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s))
                    for s in secs
                )

        resource_inputs = []
        if with_resource_overrides:
            resource_inputs = build_resource_override_maps_for_workflow(wf)
            w.inputs.extend(resource_inputs)

        # Convert self._outputs -> wdl.Output
        for o in outputs:
            sourcenode, sourcetag = o.source
            w.outputs.append(
                wdl.Output(
                    o.datatype.wdl(),
                    o.id(),
                    "{a}.{b}".format(  # Generate fully qualified stepid.tag identifier (MUST be step node)
                        a=sourcenode.id(), b=sourcetag
                    ),
                )
            )
            if o.datatype.secondary_files():
                w.outputs.extend(
                    wdl.Output(
                        o.datatype.wdl(),
                        get_secondary_tag_from_original_tag(o.id(), s),
                        "{a}.{b}".format(  # Generate fully qualified stepid.tag identifier (MUST be step node)
                            a=sourcenode.id(),
                            b=get_secondary_tag_from_original_tag(sourcetag, s),
                        ),
                    )
                    for s in o.datatype.secondary_files()
                )

        # Generate import statements (relative tool dir is later?)
        uniquetoolmap: Dict[str, Tool] = {t.versioned_id(): t for t in tools}
        w.imports = [
            wdl.Workflow.WorkflowImport(
                t.versioned_id(),
                tool_aliases[t.versioned_id().lower()].upper(),
                None if is_nested_tool else "tools/",
            )
            for t in uniquetoolmap.values()
        ]

        # Step[] -> (wdl.Task | wdl.Workflow)[]
        forbiddenidentifiers = set(
            [i.id() for i in inputs]
            + list(tool_aliases.keys())
            + list(s.id() for s in steps)
        )
        for s in steps:
            t = s.tool

            if t.versioned_id() not in wtools:
                if isinstance(t, Workflow):
                    wf_wdl, wf_tools = cls.translate_workflow(
                        t,
                        with_container=with_container,
                        is_nested_tool=True,
                        with_resource_overrides=with_resource_overrides,
                        allow_empty_container=allow_empty_container,
                        container_override=container_override,
                    )
                    wtools[t.versioned_id()] = wf_wdl
                    wtools.update(wf_tools)

                elif isinstance(t, CommandTool):
                    wtools[t.versioned_id()] = cls.translate_tool_internal(
                        t,
                        with_container=with_container,
                        with_resource_overrides=with_resource_overrides,
                        allow_empty_container=allow_empty_container,
                        container_override=container_override,
                    )
                elif isinstance(t, CodeTool):
                    wtools[t.versioned_id()] = cls.translate_code_tool_internal(
                        t,
                        with_docker=with_container,
                        with_resource_overrides=with_resource_overrides,
                        allow_empty_container=allow_empty_container,
                        container_override=container_override,
                    )

            resource_overrides = {}

            if with_resource_overrides:
                toolroverrides = build_resource_override_maps_for_tool(t)
                for r in toolroverrides:
                    resource_overrides[r.name] = s.id() + "_" + r.name

            call = translate_step_node(
                s,
                tool_aliases[t.versioned_id().lower()].upper() + "." + t.id(),
                resource_overrides,
                forbiddenidentifiers,
            )

            w.calls.append(call)

        return w, wtools

    @classmethod
    def translate_tool_inputs(cls, toolinputs: List[ToolInput]) -> List[wdl.Input]:
        ins = []
        for i in toolinputs:
            wd = i.input_type.wdl(has_default=i.default is not None)
            expr = None
            if isinstance(i.input_type, Filename):
                expr = None
            if isinstance(wd, list):
                ins.extend(wdl.Input(w, i.id()) for w in wd)
            else:

                ins.append(wdl.Input(wd, i.id(), expr, requires_quotes=False))

                sec = value_or_default(
                    i.input_type.subtype().secondary_files()
                    if isinstance(i.input_type, Array)
                    else i.input_type.secondary_files(),
                    default=[],
                )
                ins.extend(
                    wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s))
                    for s in sec
                )
        return ins

    @classmethod
    def translate_tool_outputs(
        cls, tooloutputs: List[ToolOutput], inputsmap: Dict[str, ToolInput], toolid
    ):
        outs: List[wdl.Output] = []

        for o in tooloutputs:
            outs.extend(
                translate_output_node_with_glob(o, o.glob, inputsmap, toolId=toolid)
            )
        return outs

    @classmethod
    def translate_tool_args(
        cls, toolargs: List[ToolArgument], inpmap: Dict[str, ToolInput], **debugkwargs
    ):
        if not toolargs:
            return []
        commandargs = []
        for a in toolargs:
            val = get_input_value_from_potential_selector_or_generator(a.value, inpmap)
            should_wrap_in_quotes = isinstance(val, str) and (
                a.shell_quote is None or a.shell_quote
            )
            wrapped_val = f"'{val}'" if should_wrap_in_quotes else val
            commandargs.append(
                wdl.Task.Command.CommandArgument(a.prefix, wrapped_val, a.position)
            )
        return commandargs

    @classmethod
    def build_command_from_inputs(cls, toolinputs: List[ToolInput]):
        inputsdict = {t.id(): t for t in toolinputs}
        command_ins = []
        for i in toolinputs:
            cmd = translate_command_input(i, inputsdict=inputsdict)
            if cmd:
                command_ins.append(cmd)
        return command_ins

    @classmethod
    def translate_tool_internal(
        cls,
        tool: CommandTool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override=None,
    ):

        if not Validators.validate_identifier(tool.id()):
            raise Exception(
                f"The identifier '{tool.id()}' for class '{tool.__class__.__name__}' was not validated by "
                f"'{Validators.identifier_regex}' (must start with letters, and then only contain letters, "
                f"numbers or an underscore)"
            )

        inputs: List[ToolInput] = [*cls.get_resource_override_inputs(), *tool.inputs()]
        toolouts = tool.outputs()
        inmap = {i.tag: i for i in inputs}

        if len(inputs) != len(inmap):
            dups = ", ".join(find_duplicates(list(inmap.keys())))
            raise Exception(
                f"There are {len(dups)} duplicate values in  {tool.id()}'s inputs: {dups}"
            )

        outdups = find_duplicates([o.id() for o in toolouts])
        if len(outdups) > 0:
            raise Exception(
                f"There are {len(outdups)} duplicate values in  {tool.id()}'s outputs: {outdups}"
            )

        ins: List[wdl.Input] = cls.translate_tool_inputs(inputs)
        outs: List[wdl.Output] = cls.translate_tool_outputs(toolouts, inmap, tool.id())
        command_args = cls.translate_tool_args(
            tool.arguments(), inmap, toolId=tool.id()
        )
        command_ins = cls.build_command_from_inputs(tool.inputs())

        commands = []

        env = tool.env_vars()
        if env:
            commands.extend(
                prepare_env_var_setters(env, inputsdict=inmap, toolid=tool.id())
            )

        for ti in tool.inputs():
            commands.extend(prepare_move_statements_for_input(ti))

        rbc = tool.base_command()
        bc = " ".join(rbc) if isinstance(rbc, list) else rbc

        commands.append(wdl.Task.Command(bc, command_ins, command_args))

        for ito in range(len(tool.outputs())):
            commands.extend(
                prepare_move_statements_for_output(toolouts[ito], outs[ito].expression)
            )

        r = wdl.Task.Runtime()
        if with_container:
            container = (
                WdlTranslator.get_container_override_for_tool(tool, container_override)
                or tool.container()
            )
            if container is not None:
                r.add_docker(container)
            elif not allow_empty_container:
                raise Exception(
                    f"The tool '{tool.id()}' did not have a container. Although not recommended, "
                    f"Janis can export empty docker containers with the parameter 'allow_empty_container=True "
                    f"or --allow-empty-container"
                )

        # These runtime kwargs cannot be optional, but we've enforced non-optionality when we create them
        r.kwargs["cpu"] = get_input_value_from_potential_selector_or_generator(
            CpuSelector(), inmap, string_environment=False, id="runtimestats"
        )
        r.kwargs["memory"] = '"~{select_first([runtime_memory, 4])}G"'

        if with_resource_overrides:
            ins.append(wdl.Input(wdl.WdlType.parse_type("String"), "runtime_disks"))
            r.kwargs["disks"] = "runtime_disks"
            r.kwargs["zones"] = '"australia-southeast1-b"'

        r.kwargs["preemptible"] = 2

        return wdl.Task(tool.id(), ins, outs, commands, r, version="development")

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool: CodeTool,
        with_docker=True,
        with_resource_overrides=True,
        allow_empty_container=False,
        container_override=None,
    ):
        if not Validators.validate_identifier(tool.id()):
            raise Exception(
                f"The identifier '{tool.id()}' for class '{tool.__class__.__name__}' was not validated by "
                f"'{Validators.identifier_regex}' (must start with letters, and then only contain letters, "
                f"numbers or an underscore)"
            )

        ins = [
            ToolInput(
                t.id(),
                input_type=t.intype,
                prefix=f"--{t.id()}",
                default=t.default,
                doc=t.doc,
            )
            for t in tool.tool_inputs()
        ]

        tr_ins = cls.translate_tool_inputs(cls.get_resource_override_inputs() + ins)

        outs = []
        for t in tool.tool_outputs():
            if isinstance(t.outtype, Stdout):
                outs.append(ToolOutput(t.id(), output_type=t.outtype))
                continue

            outs.append(
                ToolOutput(
                    t.id(),
                    output_type=t.outtype,
                    glob=CustomGlob(f'read_json(stdout())["{t.id()}"]'),
                )
            )

        tr_outs = cls.translate_tool_outputs(outs, {}, tool.id())

        commands = []

        scriptname = tool.script_name()

        commands.append(
            wdl.Task.Command(
                f"""
cat <<EOT >> {scriptname}
{tool.prepared_script()}
EOT"""
            )
        )

        command_ins = cls.build_command_from_inputs(ins)
        bc = tool.base_command()
        bcs = " ".join(bc) if isinstance(bc, list) else bc
        commands.append(wdl.Task.Command(bcs, command_ins, []))

        r = wdl.Task.Runtime()
        if with_docker:
            container = (
                WdlTranslator.get_container_override_for_tool(tool, container_override)
                or tool.container()
            )

            if container is not None:
                r.add_docker(container)
            elif not allow_empty_container:
                raise Exception(
                    f"The tool '{tool.id()}' did not have a container. Although not recommended, "
                    f"Janis can export empty docker containers with the parameter 'allow_empty_container=True "
                    f"or --allow-empty-container"
                )

        if with_resource_overrides:
            tr_ins.append(wdl.Input(wdl.WdlType.parse_type("String"), "runtime_disks"))
            r.kwargs["disks"] = "runtime_disks"
            r.kwargs["zones"] = '"australia-southeast1-b"'

        # These runtime kwargs cannot be optional, but we've enforced non-optionality when we create them
        # r.kwargs["cpu"] = get_input_value_from_potential_selector_or_generator(
        #     CpuSelector(), {}, string_environment=False, id="runtimestats"
        # )

        r.kwargs["memory"] = '"~{select_first([runtime_memory, 4])}G"'

        return wdl.Task(tool.id(), tr_ins, tr_outs, commands, r, version="development")

    @classmethod
    def build_inputs_file(
        cls,
        tool,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
    ) -> Dict[str, any]:
        """
        Recursive is currently unused, but eventually input overrides could be generated the whole way down
        a call chain, including subworkflows: https://github.com/openwdl/wdl/issues/217
        :param merge_resources:
        :param recursive:
        :param tool:
        :return:
        """
        from janis_core.workflow.workflow import Workflow

        inp = {}
        values_provided_from_tool = {}
        is_workflow = isinstance(tool, Workflow)

        if is_workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value or i.default
            }

        ad = {**values_provided_from_tool, **(additional_inputs or {})}

        for i in tool.tool_inputs():

            inp_key = f"{tool.id()}.{i.id()}" if is_workflow else i.id()
            value = ad.get(i.id())
            if cls.inp_can_be_skipped(i, value):
                continue

            inp_val = value

            inp[inp_key] = inp_val
            if i.intype.secondary_files():
                for sec in i.intype.secondary_files():
                    inp[
                        get_secondary_tag_from_original_tag(inp_key, sec)
                    ] = apply_secondary_file_format_to_filename(inp_val, sec)
            elif isinstance(i.intype, Array) and i.intype.subtype().secondary_files():
                # handle array of secondary files
                for sec in i.intype.subtype().secondary_files():
                    inp[get_secondary_tag_from_original_tag(inp_key, sec)] = (
                        [
                            apply_secondary_file_format_to_filename(iinp_val, sec)
                            for iinp_val in inp_val
                        ]
                        if inp_val
                        else None
                    )

        if merge_resources:
            inp.update(
                cls.build_resources_input(
                    tool, hints, max_cores, max_mem, inputs=ad, is_root=True
                )
            )

        return inp

    @classmethod
    def build_resources_input(
        cls,
        tool,
        hints,
        max_cores=None,
        max_mem=None,
        inputs=None,
        prefix=None,
        is_root=False,
    ):
        from janis_core.workflow.workflow import Workflow

        is_workflow = isinstance(tool, Workflow)
        d = super().build_resources_input(
            tool=tool,
            hints=hints,
            max_cores=max_cores,
            max_mem=max_mem,
            prefix=prefix or "",
            inputs=inputs,
        )
        if is_workflow and is_root:
            return {f"{tool.id()}.{k}": v for k, v in d.items()}
        return d

    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".wdl"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + "-inp.json"

    @staticmethod
    def tool_filename(tool):
        return (tool.versioned_id() if isinstance(tool, Tool) else str(tool)) + ".wdl"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"


def resolve_tool_input_value(tool_input: ToolInput, inputsdict, **debugkwargs):
    name = tool_input.id()
    indefault = (
        tool_input.input_type
        if isinstance(tool_input.input_type, Filename)
        else tool_input.default
    )

    default = None
    if isinstance(indefault, CpuSelector):
        if indefault.default:
            default = f"select_first([runtime_cpu, {str(indefault.default)}])"
        else:
            default = "runtime_cpu"

    elif isinstance(indefault, InputSelector):
        Logger.critical(
            f"WDL does not support command line level defaults that select a different input, this will remove the "
            f"value: '{indefault}' for tool_input '{tool_input.tag}'"
        )

    else:
        default = get_input_value_from_potential_selector_or_generator(
            indefault, inputsdict=inputsdict, string_environment=False, **debugkwargs
        )

    if default is not None:
        # Default should imply optional input
        name = f"select_first([{name}, {default}])"

    if tool_input.localise_file:
        if isinstance(tool_input.input_type, Array):
            raise Exception(
                "Localising files through `basename(x)` is unavailable for arrays of files: https://github.com/openwdl/wdl/issues/333"
            )
        name = "basename(%s)" % name

    return name


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_input.position is not None or tool_input.prefix):
        return None

    name = resolve_tool_input_value(tool_input, inputsdict=inputsdict, **debugkwargs)
    optional = tool_input.input_type.optional or (
        isinstance(tool_input.default, CpuSelector) and not tool_input.default
    )
    position = tool_input.position
    separate_value_from_prefix = tool_input.separate_value_from_prefix
    prefix = tool_input.prefix
    true = None
    sep = tool_input.separator

    is_array = isinstance(tool_input.input_type, Array)
    separate_arrays = is_array and tool_input.prefix_applies_to_all_elements

    if isinstance(tool_input.input_type, Boolean):
        true = tool_input.prefix
        prefix = None

    return wdl.Task.Command.CommandInput(
        name=name,
        optional=optional,
        prefix=prefix,
        position=position,
        separate_value_from_prefix=(
            separate_value_from_prefix
            if separate_value_from_prefix is not None
            else True
        ),
        # Instead of using default, we'll use the ~{select_first([$var, default])}
        #       (previously: ~{if defined($var) then val1 else val2})
        # as it progress through the rest properly
        # default=default,
        true=true,
        separator=(
            None
            if not is_array or separate_arrays
            else (sep if sep is not None else " ")
        ),
        separate_arrays=separate_arrays,
    )


def translate_output_node_with_glob(
    o, glob, inputmap: Dict[str, ToolInput], **debugkwargs
) -> List[wdl.Output]:
    if isinstance(o.output_type, Stdout):
        base_expression = "stdout()"
        return [wdl.Output(o.output_type.wdl(), o.id(), base_expression)]

    if isinstance(o.output_type, Stderr):
        base_expression = "stderr()"
        return [wdl.Output(o.output_type.wdl(), o.id(), base_expression)]

    elif isinstance(glob, InputSelector):
        return translate_input_selector_for_output(o, glob, inputmap, **debugkwargs)

    elif isinstance(glob, StringFormatter):
        return translate_string_formatter_for_output(o, glob, inputmap, **debugkwargs)

    elif isinstance(glob, WildcardSelector):
        base_expression = translate_wildcard_selector(glob)
        if not isinstance(o.output_type, Array):
            Logger.warn(
                f"The command tool ({debugkwargs}).{o.tag}' used a star-bind (*) glob to find the output, "
                f"but the return type was not an array. For WDL, the first element will be used, "
                f"ie: '{base_expression}[0]'"
            )
            base_expression += "[0]"
        wdl_type = wdl.WdlType.parse_type(o.output_type.wdl())
        outputs = [wdl.Output(wdl_type, o.id(), base_expression)]

        secondary = o.output_type.secondary_files()
        if secondary:
            outputs.extend(
                wdl.Output(
                    wdl_type,
                    get_secondary_tag_from_original_tag(o.id(), s),
                    base_expression,
                )
                for s in o.output_type.secondary_files()
            )
        return outputs

    elif isinstance(glob, CustomGlob):
        return [wdl.Output(o.output_type.wdl(), o.id(), glob.expression)]

    else:
        raise Exception(
            f"Tool ({debugkwargs}) has the non-selector glob: '{glob}', this is deprecated. "
            f"Please use the WildcardSelector to build output for '{o.id()}'"
        )


def translate_input_selector_for_output(
    out, selector: InputSelector, inp_map: Dict[str, ToolInput], **debugkwargs
) -> List[wdl.Output]:
    base_expression = translate_input_selector(
        selector, inp_map, string_environment=False, **debugkwargs
    )

    tool_in = inp_map.get(selector.input_to_select)
    if not tool_in:
        raise Exception(
            f"The InputSelector for tool '{debugkwargs}.{out.id()}' did not select an input (tried: '{selector.input_to_select}')"
        )
    use_basename = (
        tool_in.localise_file
        or isinstance(selector, InputSelector)
        and selector.use_basename
    )
    expression = base_expression if not use_basename else f"basename({base_expression})"

    outputs = [wdl.Output(out.output_type.wdl(), out.id(), expression)]
    for s in value_or_default(out.output_type.secondary_files(), []):
        sec_expression = None
        if "^" not in s:
            # do stuff here
            sec_expression = f'({expression}) + "{s.replace("^", "")}"'

        elif isinstance(tool_in.input_type, Filename) and tool_in.input_type.extension:
            # use the wdl function: sub
            sec_expression = 'sub({inp}, "\\\\{old_ext}$", "{new_ext}")'.format(
                inp=expression,
                old_ext=tool_in.input_type.extension,
                new_ext=s.replace("^", ""),
            )

        elif (
            File().can_receive_from(tool_in.input_type)
            and isinstance(tool_in.input_type, File)
            and tool_in.input_type.extension
        ):
            # use basename
            sec_expression = f'basename({expression}, "{tool_in.input_type.extension}") + "{s.replace("^", "")}"'

        outputs.append(
            wdl.Output(
                out.output_type.wdl(),
                get_secondary_tag_from_original_tag(out.id(), s),
                sec_expression,
            )
        )
    return outputs


def translate_string_formatter_for_output(
    out, selector: StringFormatter, inp_map: Dict[str, ToolInput], **debugkwargs
) -> List[wdl.Output]:
    """
    The output glob was a string formatter, so we'll need to build the correct glob
    by resolving the string formatter. Some minor complications involve how an output
    with a secondary file must resolve input selectors.

    For example, if you are generating an output with the Filename class, and your output
    type has secondary files, you will need to translate the generated filename with
    respect to the secondary file extension. Or the File class also has a recommended
    "extension" property now that this should consider.

    :param inp_map:
    :param out:
    :param selector:
    :return:
    """
    inputs_to_retranslate = {
        k: v
        for k, v in selector.kwargs.items()
        if not any(isinstance(v, t) for t in StringFormatter.resolved_types)
    }

    has_secondary_files = bool(out.output_type.secondary_files())

    if not has_secondary_files:
        resolved_kwargs = {
            **selector.kwargs,
            **{
                k: get_input_value_from_potential_selector_or_generator(
                    v, inputsdict=inp_map, string_environment=True, **debugkwargs
                )
                for k, v in inputs_to_retranslate.items()
            },
        }

        resolved_exp = selector.resolve_with_resolved_values(**resolved_kwargs)
        return [
            wdl.Output(
                data_type=out.output_type.wdl(),
                name=out.id(),
                expression=f'"{resolved_exp}"',
            )
        ]

    translated_inputs = {}
    input_selectors_with_secondaries = {}

    for k, v in inputs_to_retranslate.items():

        if has_secondary_files and isinstance(v, InputSelector):
            # Handle input selector separately
            input_selectors_with_secondaries[k] = v

        translated_inputs[k] = get_input_value_from_potential_selector_or_generator(
            v, inputsdict=inp_map, string_environment=True, **debugkwargs
        )

    if len(input_selectors_with_secondaries) > 1:
        invalid_keys = ", ".join(input_selectors_with_secondaries.keys())
        raise Exception(
            f"There might be an error when building the string formatter for output '{out.id()}', the "
            f"values to replace for the keys ({invalid_keys}) each which "
            f"had secondary files. This behaviour is currently undefined."
        )
    else:

        inputsel_key = next(iter(input_selectors_with_secondaries.keys()))
        inputsel = input_selectors_with_secondaries[inputsel_key]

        tool_in = inp_map.get(inputsel.input_to_select)

    expression = selector.resolve_with_resolved_values(
        **{**selector.kwargs, **translated_inputs}
    )

    outputs = [wdl.Output(out.output_type.wdl(), out.id(), f'"{expression}"')]
    for s in value_or_default(out.output_type.secondary_files(), []):
        sec_expression = None
        if "^" not in s:
            # do stuff here
            sec_expression = expression + s

        elif tool_in:
            if isinstance(tool_in.input_type, Filename):
                if not tool_in.input_type.extension:
                    raise Exception(
                        f"Unsure how to handle secondary file '{s}' as it uses the escape characater '^' but"
                        f"Janis can't determine the extension of the input '{tool_in.id()} to replace in the "
                        f"WDL translation. You will to annotate the input with an extension Filename(extension=)"
                    )

                # use the wdl function: sub
                sec_expression = '~{sub({inp}, "\\\\{old_ext}$", "{new_ext}")}'.format(
                    inp=expression,
                    old_ext=tool_in.input_type.extension,
                    new_ext=s.replace("^", ""),
                )

            elif File().can_receive_from(tool_in.input_type):
                if not tool_in.input_type.extension:
                    raise Exception(
                        f"Unsure how to handle secondary file '{s}' as it uses the escape characater '^' but"
                        f"Janis can't determine the extension of the input '{tool_in.id()} to replace in the "
                        f"WDL translation. When creating an instance of the {tool_in.input_type.__name__} class, "
                        f"it's possible to annotate an extension. For example, the Zip filetype is going to "
                        f'have the .zip extension, so could be initialised with Zip(extension=".zip"'
                    )
                # use basename
                replaced_s = s.replace("^", "")
                sec_expression = f'basename({tool_in.id()}, "{tool_in.input_type.extension}") + "{replaced_s}"'
            else:
                raise Exception(
                    f"Unsure how to handle secondary file '{s}' as it uses the escape characater '^' but"
                    f"Janis can't determine the extension of the input '{tool_in.id()} to replace in the "
                    f"WDL translation. Extra handling in 'wdl.translate_string_formatter_for_output' may be"
                    f"required for the '{tool_in.input_type.__name__}' type might be required."
                )
        else:
            sec_expression = apply_secondary_file_format_to_filename(expression, s)

        outputs.append(
            wdl.Output(
                out.output_type.wdl(),
                get_secondary_tag_from_original_tag(out.id(), s),
                f'"{sec_expression}"',
            )
        )

    return outputs


def translate_step_node(
    node2,
    step_identifier: str,
    resource_overrides: Dict[str, str],
    invalid_identifiers: Set[str],
) -> wdl.WorkflowCallBase:
    """
    Convert a step into a wdl's workflow: call { **input_map }, this handles creating the input map and will
    be able to handle multiple scatters on this step node. If there are multiple scatters, the scatters will be ordered
    in to out by alphabetical order.
    
    This method isn't perfect, when there are multiple sources it's not correctly resolving defaults,
    and tbh it's pretty confusing.
    
    :param node:
    :param step_identifier:
    :param step_alias:
    :param resource_overrides:
    :return:
    """
    from janis_core.workflow.workflow import StepNode, InputNode

    node: StepNode = node2
    step_alias: str = node.id()

    ins = node.inputs()

    # Sanity check our step node connections:

    # 1. make sure our inputs are all present:
    missing_keys = [
        k
        for k in ins.keys()
        if k not in node.sources and not (ins[k].intype.optional or ins[k].default)
    ]
    if missing_keys:
        raise Exception(
            f"Error when building connections for step '{node.id()}', "
            f"missing the required connection(s): '{', '.join(missing_keys)}'"
        )

    # 2. gather the scatters, and make sure none of them are derived from multiple sources, otherwise
    #       we're like double zipping things, it's complicated and it's even MORE complicated in WDL.
    scatterable: List[StepTagInput] = []

    if node.scatter:
        unbound_scatter_keys = [k for k in node.scatter.fields if k not in node.sources]
        if len(unbound_scatter_keys):
            raise Exception(
                f"Attempted to scatter {node.id()} on field(s) [{', '.join(unbound_scatter_keys)}] however "
                "these inputs were not mapped on step construction. Make sure that those unbound keys exist"
                f"in your step definition (eg: "
                f"{node.tool.__class__.__name__}({', '.join(k + '=inp' for k in unbound_scatter_keys)})"
            )
        scatterable = [node.sources[k] for k in node.scatter.fields]

        invalid_sources = [
            si
            for si in scatterable
            if si.multiple_inputs or isinstance(si.source(), list)
        ]
        if len(invalid_sources) > 0:
            invalid_sources_str = ", ".join(
                si.dotted_source() for si in invalid_sources
            )
            raise NotImplementedError(
                f"The edge(s) '{invalid_sources_str}' on node '{node.id()}' scatters"
                f"on multiple inputs, this behaviour has not been implemented"
            )

    # 1. Generate replacement of the scatterable key(s) with some random variable, eg: for i in iterable:
    scattered_old_to_new_identifier = generate_scatterable_details(
        scatterable, forbiddenidentifiers=invalid_identifiers
    )

    # Let's map the inputs, to the source. We're using a dictionary for the map atm, but WDL requires the _format:
    #       fieldName: sourceCall.Output

    inputs_map = {}
    for k in ins:
        if k not in node.sources:
            continue

        array_input_from_single_source = False

        edge: StepTagInput = node.sources[k]
        source: Edge = edge.source()  # potentially single item or array

        # We have multiple sources going to the same entry
        if isinstance(source, list):
            if len(source) == 1:
                source = source[0]
            elif len(source) > 1:

                # There are multiple sources, this is sometimes a little tricky
                input_name_maps = ", ".join(edge.dotted_source())

                multiple_sources_failure_reasons = []

                unique_types = set()
                for x in edge.source():
                    t: TOutput = (
                        first_value(x.start.outputs())
                        if not x.stag
                        else x.start.outputs()[x.stag]
                    )

                    unique_types.update(t.outtype.secondary_files() or [""])
                if len(unique_types) > 1:
                    multiple_sources_failure_reasons.append(
                        f"has {len(unique_types)} different DataTypes with varying secondaries"
                    )
                if node.scatter:
                    multiple_sources_failure_reasons.append(f"is scattered")

                if len(multiple_sources_failure_reasons) > 0:
                    reasons = " and ".join(multiple_sources_failure_reasons)
                    Logger.critical(
                        f"Conversion to WDL for field '{node.id()}.{k}' does not fully support multiple sources."
                        f" This will only work if all of the inputs ({input_name_maps}) have the same secondaries "
                        f"AND this field ('{k}') is not scattered. However this connection {reasons}"
                    )

                ds = edge.dotted_source()
                if isinstance(ds, list):
                    ds = "[" + ", ".join(ds) + "]"
                inputs_map[k] = ds
                f = edge.finish.inputs()[edge.ftag]
                secs = (
                    f.intype.subtype().secondary_files()
                    if isinstance(f.intype, Array)
                    else f.intype.secondary_files()
                )
                if secs:
                    for sec in secs:
                        inputs_map[get_secondary_tag_from_original_tag(k, sec)] = (
                            "["
                            + ", ".join(
                                get_secondary_tag_from_original_tag(kk, sec)
                                for kk in edge.dotted_source()
                            )
                            + "]"
                        )
                continue

        elif source:
            it = source.finish.inputs()[source.ftag].intype
            if source.stag:
                ot = source.start.outputs()[source.stag].outtype
            else:
                ot = first_value(source.start.outputs()).outtype
            if (
                isinstance(it, Array)
                and not isinstance(ot, Array)
                and not source.scatter
            ):
                array_input_from_single_source = True
        secondary = None
        # We're connecting to another step
        if source and isinstance(source.finish, StepNode) and source.ftag:

            it = source.finish.inputs()[source.ftag].intype

            secondary = (
                it.subtype().secondary_files()
                if isinstance(it, Array)
                else it.secondary_files()
            )
            if secondary and isinstance(source.start, StepNode) and source.stag:

                ot = source.start.outputs()[source.stag].outtype

                sec_out = set(
                    value_or_default(
                        ot.subtype().secondary_files()
                        if isinstance(ot, Array)
                        else ot.secondary_files(),
                        default=[],
                    )
                )
                sec_in = set(secondary)
                if not sec_out.issubset(sec_in):
                    raise Exception(
                        f"An error occurred when connecting '{source.source_dotted()}' to "
                        f"'{source.finish.id()}.{source.ftag}', there were secondary files in the final node "
                        f"that weren't present in the source: {', '.join(sec_out.difference(sec_in))}"
                    )

        # Edge defaults
        if not source:
            # edge but no source, probably a default
            if not edge.default:
                Logger.critical(
                    f"Skipping connection to '{edge.finish}.{edge.ftag}' had no source or default, "
                    f"please raise an issue as investigation may be required"
                )
                continue

            if isinstance(edge.default, bool):
                inputs_map[k] = "true" if edge.default else "false"
            elif isinstance(edge.default, str):
                inputs_map[k] = f'"{edge.default}"'
            else:
                inputs_map[k] = edge.default

        # Scattering on multiple secondary files
        elif edge in scatterable and secondary:
            # We're ensured through inheritance and .receiveBy that secondary files will match.
            ds = source.source_dotted()
            Logger.log(
                f"Oh boii, we're gonna have some complicated scattering here with {len(secondary)} secondary file(s)"
            )

            identifier = scattered_old_to_new_identifier[ds]
            inputs_map[k] = identifier[0] + "[0]"
            for idx in range(len(secondary)):
                sec = secondary[idx]
                inputs_map[
                    get_secondary_tag_from_original_tag(k, sec)
                ] = f"{identifier[0]}[{idx + 1}]"

        else:
            ds = source.source_dotted()
            default = None
            if source.start and isinstance(source.start, InputNode):
                default = source.start.default

            inpsourcevalue = None
            if (
                ds in scattered_old_to_new_identifier
                and scattered_old_to_new_identifier[ds]
            ):
                # can't get here with secondary

                s = scattered_old_to_new_identifier[ds]
                inpsourcevalue = s[0]
                default = None

            else:
                inpsourcevalue = ds
                if secondary:
                    if default:
                        Logger.critical(
                            f"The default '{default}' will not be applied to the map '{ds}' as there are secondary files"
                        )
                        default = None
                    for idx in range(len(secondary)):
                        sec = secondary[idx]
                        inputs_map[
                            get_secondary_tag_from_original_tag(k, sec)
                        ] = get_secondary_tag_from_original_tag(ds, sec)

            if default:
                defval = get_input_value_from_potential_selector_or_generator(
                    default, inputsdict=None, string_environment=False, dottedsource=ds
                )

                if isinstance(defval, bool):
                    defval = "true" if defval else "false"

                inpsourcevalue = f"select_first([{inpsourcevalue}, {defval}])"

            if array_input_from_single_source and not (
                isinstance(source.start, StepNode) and source.start.scatter
            ):
                inpsourcevalue = f"[{inpsourcevalue}]"
                if secondary:
                    for sec in secondary:
                        tag = get_secondary_tag_from_original_tag(k, sec)
                        inputs_map[tag] = f"[{inputs_map[tag]}]"

            inputs_map[k] = inpsourcevalue

    inputs_map.update(resource_overrides)

    call = wdl.WorkflowCall(step_identifier, step_alias, inputs_map)

    if len(scatterable) > 0:
        call = wrap_scatter_call(
            call, node.scatter, scatterable, scattered_old_to_new_identifier
        )

    return call


def generate_scatterable_details(
    scatterable: List[StepTagInput], forbiddenidentifiers: Set[str]
):

    # this dictionary is what we're going to use to map our current
    # identifier to the scattered identifier. This step is just the
    # setup, and in the next for loop, we'll
    scattered_old_to_new_identifier = {
        k.dotted_source(): (k.dotted_source(), k.source())
        for k in scatterable
        if not isinstance(k.dotted_source(), list)
    }

    # Make a copy of the forbiddenIds and add the identifiers of the source
    forbiddenidentifierscopy = set(forbiddenidentifiers).union(
        set(v[0] for v in scattered_old_to_new_identifier.values())
    )

    # We'll wrap everything in the scatter block later, but let's replace the fields we need to scatter
    # with the new scatter variable (we'll try to guess one based on the fieldname). We might need to eventually
    # pass the workflow inputs to make sure now conflict will arise.

    if len(scatterable) > 1:
        # We'll generate one variable in place
        standin = generate_new_id_from("Q", forbiddenidentifierscopy)

        # Store the the standin variable for us to use later (as an illegal identifier character)
        scattered_old_to_new_identifier["-"] = standin
        forbiddenidentifierscopy.add(standin)

        # Then we'll need to generate something like:
        #       A -> Q[0], B -> Q[0][0] -> ..., N -> Q([0] * (n-1))[1]
        n = len(scatterable)
        for i in range(len(scatterable)):
            s = scatterable[i]
            newid = standin + (i) * ".right" + ".left"
            if i == n - 1:
                newid = standin + (n - 1) * ".right"
            scattered_old_to_new_identifier[s.dotted_source()] = (
                newid,
                first_value(s.source_map).start,
            )
            forbiddenidentifierscopy.add(newid)
    else:

        for s in scatterable:

            # We asserted earlier that the source_map only has one value (through multipleInputs)
            e: Edge = first_value(s.source_map)
            newid = generate_new_id_from(
                e.stag or e.source_dotted().replace(".", ""), forbiddenidentifierscopy
            )
            scattered_old_to_new_identifier[s.dotted_source()] = (newid, e.start)
            forbiddenidentifierscopy.add(newid)

    return scattered_old_to_new_identifier


def wrap_scatter_call(
    call, scatter: ScatterDescription, scatterable, scattered_old_to_new_identifier
):
    from janis_core.workflow.workflow import InputNode

    # Let's start the difficult process of scattering, in WDL we'll:
    #
    #       1. We already generated the new "scatter" variable that will
    #            be used in place of the original .dotted_source()
    #
    #       2. Generate annotations for accessory / secondary files
    #
    #       3. Wrap everything in the appropriate scatter block,
    #           especially considering scattering by multiple variables
    #

    # 2. Explanation:
    #     So, the current way of mixing accessory files is not really                        _
    #     supported, but a little complicated basically, if our scatterable edge           _| |
    #     contains secondary files, they'll all be arrays of separate files, eg:         _| | |
    #                                                                                   | | | | __
    #     File[] bams =     [...]                                                       | | | |/  \
    #     File[] bams_bai = [...]                                                       |       /\ \
    #                                                                                   |       \/ /
    #     We can handle this by transposing the array of both items, eg:                 \        /
    #                                                                                     |      /
    #         transpose([bam1, bam2, ..., bamn], [bai1, bai2, ..., bai3])                 |     |
    #               => [[bam1, bai1], [bam2, bai2], ..., [bamn, bain]]
    #
    #     and then unwrap them using their indices and hopefully have everything line up:
    #
    #     Source: https://software.broadinstitute.org/wdl/documentation/spec#arrayarrayx-transposearrayarrayx

    # We need to generate the source statement, this is easy if we're scattering by one variable

    # sanity check
    if len(scatterable) == 0:
        return call

    # generate the new source map

    insource_ar = []
    for s in scatterable:
        secondary = s.finish.tool.inputs_map()[s.ftag].intype.secondary_files()
        if secondary:
            ds = s.dotted_source()
            joined_tags = ", ".join(
                get_secondary_tag_from_original_tag(ds, sec) for sec in secondary
            )
            transformed = f"transpose([{ds}, {joined_tags}])"
            insource_ar.append(transformed)

        else:
            (newid, startnode) = scattered_old_to_new_identifier[s.dotted_source()]
            insource = s.dotted_source()
            if isinstance(startnode, InputNode) and startnode.default is not None:
                resolved = get_input_value_from_potential_selector_or_generator(
                    startnode.default,
                    None,
                    string_environment=False,
                    scatterstep=insource,
                )
                if isinstance(resolved, bool):
                    resolved = "true" if resolved else "false"

                insource_ar.append(f"select_first([{insource}, {resolved}])")
            else:
                insource_ar.append(insource)

    insource = None
    alias = None
    if len(insource_ar) == 1:
        insource = insource_ar[0]
        alias = first_value(scattered_old_to_new_identifier)[0]
    else:
        method = "zip" if scatter.method == ScatterMethods.dot else "cross"
        insource = recursive_2param_wrap(method, insource_ar)
        alias = scattered_old_to_new_identifier["-"]

    return wdl.WorkflowScatter(alias, insource, [call])


## SELECTOR HELPERS


def get_input_value_from_potential_selector_or_generator(
    value, inputsdict, string_environment=True, **debugkwargs
):
    """
    We have a value which could be anything, and we want to convert it to a expressionable value.
    If we have a string, it should be in quotes, etc. It should be "Type paramname = <expressionable>"
    :param value:
    :param tool_id:
    :param string_environment:  Do we need to wrap string literals in quotes,
                                or do we need to wrap variables in expr block.
    :return:
    """
    if value is None:
        return None

    if isinstance(value, list):
        toolid = value_or_default(debugkwargs.get("tool_id"), "get-value-list")
        joined_values = ", ".join(
            str(
                get_input_value_from_potential_selector_or_generator(
                    value[i],
                    inputsdict,
                    string_environment=False,
                    tool_id=toolid + "." + str(i),
                )
            )
            for i in range(len(value))
        )
        return f"[{joined_values}]"
    elif isinstance(value, str):
        return value if string_environment else f'"{value}"'
    elif isinstance(value, bool):
        return "true" if value else "false"
    elif isinstance(value, int) or isinstance(value, float):
        return value
    elif isinstance(value, Filename):
        gen_filename = value.generated_filename(
            inputs=prepare_filename_replacements_for(value.prefix, inputsdict)
        )
        return gen_filename if string_environment else f'"{gen_filename}"'
    elif isinstance(value, StringFormatter):
        return translate_string_formatter(
            selector=value,
            inputsdict=inputsdict,
            string_environment=string_environment,
            **debugkwargs,
        )
    elif isinstance(value, WildcardSelector):
        raise Exception(
            f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
        )
    # elif isinstance(value, CpuSelector):
    #     return translate_cpu_selector(
    #         value, inputsdict, string_environment=string_environment
    #     )
    # elif isinstance(value, MemorySelector):
    #     return translate_mem_selector(
    #         value, inputsdict, string_environment=string_environment
    #     )
    elif isinstance(value, InputSelector):
        return translate_input_selector(
            selector=value,
            inputsdict=inputsdict,
            string_environment=string_environment,
            **debugkwargs,
        )
    elif callable(getattr(value, "wdl", None)):
        return value.wdl()

    warning = ""
    if isclass(value):
        stype = value.__name__
        warning = f", this is likely due to the '{stype}' not being initialised"
    else:
        stype = value.__class__.__name__
    raise Exception(
        f"Could not detect type '{stype}' to convert to input value{warning}"
    )


def translate_string_formatter(
    selector: StringFormatter, inputsdict, string_environment=False, **debugkwargs
):
    # we should raise an Exception if any of our inputs are optional without a default

    invalid_select_inputs = [
        (k, selector.kwargs[k].input_to_select)
        for k in selector.kwargs
        # Our selector is getting an input
        if isinstance(selector.kwargs[k], InputSelector)
        and selector.kwargs[k].input_to_select in inputsdict
        and not isinstance(
            inputsdict[selector.kwargs[k].input_to_select].input_type, Filename
        )
        # our selected input is optional
        and inputsdict[selector.kwargs[k].input_to_select].input_type.optional
        # our selected input does NOT have a default
        # tbh, this ToolInput might have a default that selects a different input that is null,
        # but I'm not going down this rabbit hole
        and inputsdict[selector.kwargs[k].input_to_select].default is None
    ]

    if len(invalid_select_inputs) > 0:
        tags = ", ".join(f"'{k[0]}'" for k in invalid_select_inputs)
        inps = ", ".join(f"'{k[1]}'" for k in invalid_select_inputs)
        raise Exception(
            f'There was an error when resolving the format "{selector._format}", the tag(s) {tags} respectively '
            f"selected input(s) {inps} that were optional and did NOT have a default value."
        )

    value = selector.resolve_with_resolved_values(
        **{
            k: get_input_value_from_potential_selector_or_generator(
                selector.kwargs[k],
                inputsdict=inputsdict,
                string_environment=True,
                **debugkwargs,
            )
            for k in selector.kwargs
        }
    )

    return value if string_environment else f'"{value}"'


def translate_input_selector(
    selector: InputSelector, inputsdict, string_environment=True, **debugkwargs
):
    if not selector.input_to_select:
        raise Exception("No input was selected for input selector: " + str(selector))

    if not inputsdict:
        raise Exception(
            f"An internal error has occurred when selecting the input '{selector.input_to_select}'"
        )

    if selector.input_to_select not in inputsdict:
        raise Exception(
            f"Couldn't find input '{selector.input_to_select}' in tool '{debugkwargs}'"
        )

    inp = inputsdict[selector.input_to_select]
    name = resolve_tool_input_value(inp, inputsdict, **debugkwargs)
    if string_environment:
        return f"~{{{name}}}"
    else:
        return name


# These are now handled entirely by the translate_input_selector_class

# def translate_cpu_selector(selector: CpuSelector, inputsdict, string_environment=True):
#     return translate_input_selector(selector, inputsdict, string_environment=string_environment)
# def translate_mem_selector(selector: MemorySelector, inputsdict, string_environment=True):
#     return translate_input_selector(selector, inputsdict, string_environment=string_environment)


def translate_wildcard_selector(selector: WildcardSelector):
    if not selector.wildcard:
        raise Exception(
            "No wildcard was selected for wildcard selector: " + str(selector)
        )
    return f'glob("{selector.wildcard}")'


## HELPER METHODS


def value_or_default(ar, default):
    """
    default is ar is None, else return the value of ar, returns ar even if ar is FALSEY
    :param ar:
    :param default:
    :return:
    """
    return default if ar is None else ar


def build_aliases(steps2):
    """
    From a list of stepNodes, generate the toolname alias (tool name to unique import alias)
    and the step alias, which
    :param steps: list of step nodes
    :return:
    """
    from janis_core.workflow.workflow import StepNode

    steps: List[StepNode] = steps2

    get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])
    aliases: Set[str] = set()

    tools: List[Tool] = [s.tool for s in steps]
    tool_name_to_tool: Dict[str, Tool] = {t.versioned_id().lower(): t for t in tools}
    tool_name_to_alias = {}
    steps_to_alias: Dict[str, str] = {
        s.id().lower(): get_alias(s.id()).lower() for s in steps
    }

    for tool in tool_name_to_tool:
        a = get_alias(tool).upper()
        s = a
        idx = 2
        while s in aliases:
            s = a + str(idx)
            idx += 1
        aliases.add(s)
        tool_name_to_alias[tool] = s

    return tool_name_to_alias, steps_to_alias


def get_secondary_tag_from_original_tag(original, secondary):
    secondary_without_punctuation = secondary.replace(".", "").replace("^", "")
    return original + "_" + secondary_without_punctuation


def prepare_env_var_setters(
    reqs: Dict[str, Any], inputsdict, **debugkwargs
) -> List[wdl.Task.Command]:
    if not reqs:
        return []

    statements = []
    for k, v in reqs.items():
        val = get_input_value_from_potential_selector_or_generator(
            v, inputsdict=inputsdict, string_environment=True, **debugkwargs
        )
        statements.append(wdl.Task.Command(f"export {k}='{val}'"))

    return statements


def prepare_move_statements_for_input(ti: ToolInput):
    """
    Update 2019-12-16:

        ToolInput introduces 'presents_as' and 'secondaries_present_as' fields (in addition to 'localise').
        This allows you to present a file as a specific filename, or change how the secondary files present to the tool.

        Additional considerations:
            - we MUST ensure that the secondary file is in the same directory as the base
            - if something is unavailable for an array, we must let the user know

        This is the logic that should be applied:

            - localise=True :: Moves the file into the execution directory ('.')
            - presents_as :: rewrites the input file to be the requires name (in the same directory).
                             This should also rewrite the secondary files per the original rules
            - secondaries_present_as :: rewrites the extension of the secondary files

        Combinations of these can be used

            - presents_as + localise :: should move the file into the execution directory as the new name
            - secondaries_present_as + presents as :: should rewrite the secondary

        The easiest way to do this is:

            - if localise or presents_as is given, generate the new filename
            - apply the secondary file naming rules to all the others (ensures it's in the same directory for free).

    """

    it = ti.input_type
    commands: List[wdl.Task.Command] = []

    if not (ti.localise_file or ti.presents_as or ti.secondaries_present_as):
        return commands

    if not issubclass(type(it), File):
        Logger.critical(
            "Janis has temporarily removed support for localising array types"
        )
        return commands

    base = f"~{{{ti.id()}}}"

    if ti.localise_file or ti.presents_as:
        newlocation = None

        if ti.localise_file and not ti.presents_as:
            newlocation = "."
        elif not ti.localise_file and ti.presents_as:
            newlocation = f'"`dirname ~{{{ti.id()}}}`/{ti.presents_as}"'
            base = newlocation
        else:
            newlocation = ti.presents_as
            base = f'"{ti.presents_as}"'

        commands.append(wdl.Task.Command(f"ln -f ~{{{ti.id()}}} {newlocation}"))

    if it.secondary_files():
        for s in it.secondary_files():
            sectag = get_secondary_tag_from_original_tag(ti.id(), s)

            newext, iters = split_secondary_file_carats(
                ti.secondaries_present_as.get(s, s)
            )
            newpath = REMOVE_EXTENSION(base, iters) + newext
            commands.append(wdl.Task.Command(f"ln -f ~{{{sectag}}} {newpath}"))

    return commands


def prepare_move_statements_for_output(
    to: ToolOutput, baseexpression
) -> List[wdl.Task.Command]:
    """
    Update 2019-12-16:

            - presents_as :: rewrites the output file to be the requires name (in the same directory).
                             This should also rewrite the secondary files per the original rules
            - secondaries_present_as :: rewrites the extension of the secondary files

        Combinations of these can be used

            - presents_as + localise :: should move the file into the execution directory as the new name
            - secondaries_present_as + presents as :: should rewrite the secondary
    """

    ot = to.output_type
    commands = []

    if not (to.presents_as or to.secondaries_present_as):
        return commands

    if not issubclass(type(ot), File):
        Logger.critical(
            f"Janis has temporarily removed support for localising '{type(ot)}' types"
        )
        return commands

    base = f"~{{{baseexpression}}}"

    if to.presents_as:
        newlocation = to.presents_as
        base = f'"{to.presents_as}"'

        commands.append(wdl.Task.Command(f"ln -f ~{{{to.id()}}} {newlocation}"))

    if to.secondaries_present_as and ot.secondary_files():
        for s in ot.secondary_files():
            if s not in to.secondaries_present_as:
                continue

            newextvalues = split_secondary_file_carats(s)
            oldextvalues = split_secondary_file_carats(to.secondaries_present_as[s])

            oldpath = REMOVE_EXTENSION(base, oldextvalues[1]) + oldextvalues[0]
            newpath = REMOVE_EXTENSION(base, newextvalues[1]) + newextvalues[0]

            commands.append(wdl.Task.Command(f"ln -f {oldpath} {newpath}"))

    return commands


def build_resource_override_maps_for_tool(tool, prefix=None) -> List[wdl.Input]:
    inputs = []

    if not prefix:
        prefix = ""  # wf.id() + "."
    else:
        prefix += "_"

    if isinstance(tool, (CommandTool, CodeTool)):
        inputs.extend(
            [
                wdl.Input(wdl.WdlType.parse_type("Int?"), prefix + "runtime_memory"),
                wdl.Input(wdl.WdlType.parse_type("Int?"), prefix + "runtime_cpu"),
                wdl.Input(wdl.WdlType.parse_type("String"), prefix + "runtime_disks"),
            ]
        )
    else:
        inputs.extend(build_resource_override_maps_for_workflow(tool, prefix=prefix))

    return inputs


def build_resource_override_maps_for_workflow(wf, prefix=None) -> List[wdl.Input]:

    # returns a list of key, value pairs
    inputs = []

    for s in wf.step_nodes.values():
        tool: Tool = s.tool

        tool_pre = (prefix or "") + s.id()
        inputs.extend(build_resource_override_maps_for_tool(tool, prefix=tool_pre))

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
            base = f'basename({tinp.id()}, "{intype.extension}")'
        else:
            base = f"basename({tinp.id()})"
    else:
        base = tinp.id()

    if intype.optional:
        replacement = f'~{{if defined({tinp.id()}) then {base} else "generated"}}'
    else:
        replacement = f"~{{{base}}}"

    return {inp.input_to_select: replacement}
