"""
WDL

This is one of the more complicated classes, it takes the janis in-memory representation of a workflow,
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
from typing import List, Dict, Optional, Any, Set, Tuple

import wdlgen as wdl

from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.tool.commandtool import CommandTool
from janis_core.tool.tool import Tool, ToolInput, ToolArgument, ToolOutput
from janis_core.translations.translationbase import TranslatorBase
from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    MemorySelector,
    StringFormatter,
)
from janis_core.types.common_data_types import (
    Stdout,
    Array,
    Boolean,
    Filename,
    File,
    Int,
)
from janis_core.utils import first_value
from janis_core.utils.logger import Logger
from janis_core.utils.validators import Validators

# from janis_core.workflow.step import StepNode


## PRIMARY TRANSLATION METHODS


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
        cls, wfi, with_docker=True, with_resource_overrides=False, is_nested_tool=False
    ) -> Tuple[any, Dict[str, any]]:
        """
        Translate the workflow into wdlgen classes!


        :param with_resource_overrides:
        :param with_docker:
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
            wd = i.datatype.wdl(has_default=i.default is not None)

            expr = None
            if isinstance(i.datatype, Filename):
                expr = f'"{i.datatype.generated_filename()}"'

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
        uniquetoolmap: Dict[str, Tool] = {t.id(): t for t in tools}
        w.imports = [
            wdl.Workflow.WorkflowImport(
                t.id(),
                tool_aliases[t.id().lower()].upper(),
                None if is_nested_tool else "tools/",
            )
            for t in uniquetoolmap.values()
        ]

        # Step[] -> (wdl.Task | wdl.Workflow)[]
        for s in steps:
            t = s.tool

            if t.id() not in wtools:
                if isinstance(t, Workflow):
                    wf_wdl, wf_tools = cls.translate_workflow(
                        t,
                        with_docker=with_docker,
                        is_nested_tool=True,
                        with_resource_overrides=with_resource_overrides,
                    )
                    wtools[t.id()] = wf_wdl
                    wtools.update(wf_tools)

                elif isinstance(t, CommandTool):
                    wtools[t.id()] = cls.translate_tool_internal(
                        t,
                        with_docker=with_docker,
                        with_resource_overrides=with_resource_overrides,
                    )

            resource_overrides = {}
            for r in resource_inputs:
                if not r.name.startswith(s.id()):
                    continue

                resource_overrides[r.name[(len(s.id()) + 1) :]] = r.name
            call = translate_step_node(
                s,
                tool_aliases[t.id().lower()].upper() + "." + t.id(),
                s.id(),
                resource_overrides,
            )

            w.calls.append(call)

        return w, wtools

    @classmethod
    def translate_tool_inputs(cls, toolinputs: List[ToolInput]):
        ins = []
        for i in toolinputs:
            wd = i.input_type.wdl(has_default=i.default is not None)
            expr = None
            if isinstance(i.input_type, Filename):
                expr = f'"{i.input_type.generated_filename()}"'
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
    def build_command_from_inputs(
        cls, toolinputs: List[ToolInput], inputmap: Dict[str, ToolInput]
    ):
        command_ins = []
        for i in toolinputs:
            cmd = translate_command_input(i, inputmap)
            if cmd:
                command_ins.append(cmd)
        return command_ins

    @classmethod
    def translate_tool_internal(
        cls, tool: CommandTool, with_docker=True, with_resource_overrides=False
    ):

        if not Validators.validate_identifier(tool.id()):
            raise Exception(
                f"The identifier '{tool.id()}' for class '{tool.__class__.__name__}' was not validated by "
                f"'{Validators.identifier_regex}' (must start with letters, and then only contain letters, "
                f"numbers or an underscore)"
            )

        inputs: List[ToolInput] = [*cls.get_resource_override_inputs(), *tool.inputs()]
        inmap = {i.tag: i for i in inputs}

        ins: List[wdl.Input] = cls.translate_tool_inputs(inputs)
        outs: List[wdl.Output] = cls.translate_tool_outputs(
            tool.outputs(), inmap, tool.id()
        )
        command_args = cls.translate_tool_args(
            tool.arguments(), inmap, toolId=tool.id()
        )
        command_ins = cls.build_command_from_inputs(tool.inputs(), inmap)

        commands = []
        for ti in tool.inputs():
            if ti.localise_file:
                commands.extend(prepare_move_statement_for_input_to_localise(ti))
            else:
                commands.extend(prepare_move_statement_for_potential_secondaries(ti))

        rbc = tool.base_command()
        bc = " ".join(rbc) if isinstance(rbc, list) else rbc

        commands.append(wdl.Task.Command(bc, command_ins, command_args))

        r = wdl.Task.Runtime()
        if with_docker:
            r.add_docker(tool.container())

        # These runtime kwargs cannot be optional, but we've enforced non-optionality when we create them
        r.kwargs["cpu"] = get_input_value_from_potential_selector_or_generator(
            CpuSelector(), inmap, string_environment=False, id="runtimestats"
        )
        r.kwargs["memory"] = wdl.IfThenElse(
            "defined(runtime_memory)", '"${runtime_memory}G"', '"4G"'
        )

        if with_resource_overrides:
            ins.append(wdl.Input(wdl.WdlType.parse_type("String"), "runtime_disks"))
            r.kwargs["disks"] = "runtime_disks"
            r.kwargs["zones"] = '"australia-southeast1-b"'

        r.kwargs["preemptible"] = 2

        return wdl.Task(tool.id(), ins, outs, commands, r, version="development")

    @classmethod
    def build_inputs_file(
        cls,
        workflow,
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
        :param workflow:
        :return:
        """
        inp = {}
        ad = additional_inputs or {}

        for i in workflow.input_nodes.values():
            inp_key = f"{workflow.id()}.{i.id()}"
            value = ad.get(i.id()) or i.value or i.default
            if cls.inp_can_be_skipped(i, value):
                continue

            inp_val = value

            inp[inp_key] = inp_val
            if i.datatype.secondary_files():
                for sec in i.datatype.secondary_files():
                    inp[
                        get_secondary_tag_from_original_tag(inp_key, sec)
                    ] = apply_secondary_file_format_to_filename(inp_val, sec)
            elif (
                isinstance(i.datatype, Array) and i.datatype.subtype().secondary_files()
            ):
                # handle array of secondary files
                for sec in i.datatype.subtype().secondary_files():
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
                    workflow, hints, max_cores, max_mem, inputs=ad
                )
            )

        return inp

    @classmethod
    def build_resources_input(
        cls, workflow, hints, max_cores=None, max_mem=None, inputs=None, prefix=None
    ):
        return super().build_resources_input(
            workflow=workflow,
            hints=hints,
            max_cores=max_cores,
            max_mem=max_mem,
            prefix=prefix or f"{workflow.id()}.",
            inputs=inputs,
        )

    @staticmethod
    def workflow_filename(workflow):
        return workflow.id() + ".wdl"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + "-inp.json"

    @staticmethod
    def tool_filename(tool):
        return (tool.id() if isinstance(tool, Tool) else str(tool)) + ".wdl"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"


def resolve_tool_input_value(tool_input: ToolInput, **debugkwargs):
    name = tool_input.id()
    indefault = (
        tool_input.input_type
        if isinstance(tool_input.input_type, Filename)
        else tool_input.default
    )

    default = None
    if isinstance(indefault, CpuSelector):
        if indefault.default:
            default = "if defined(runtime_cpu) then runtime_cpu else " + str(
                indefault.default
            )
        else:
            default = "runtime_cpu"

    elif isinstance(indefault, InputSelector):
        Logger.critical(
            f"WDL does not support command line level defaults that select a different input, this will remove the "
            f"value: '{indefault}' for tool_input '{tool_input.tag}'"
        )

    else:
        default = get_input_value_from_potential_selector_or_generator(
            indefault, None, string_environment=False, **debugkwargs
        )

    if default:
        name = f"if defined({name}) then {name} else {default}"

    if tool_input.localise_file:
        if isinstance(tool_input.input_type, Array):
            raise Exception(
                "Localising files through `basename(x)` is unavailable for arrays of files: https://github.com/openwdl/wdl/issues/333"
            )
        name = "basename(%s)" % name

    return name


def translate_command_input(tool_input: ToolInput, inputsdict, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    # TODO: make a property on ToolInput (.bind_to_commandline) and set default to true
    if not (tool_input.position is not None or tool_input.prefix):
        return None

    name = resolve_tool_input_value(tool_input, **debugkwargs)
    optional = tool_input.input_type.optional or isinstance(
        tool_input.default, CpuSelector
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
        # Instead of using default, we'll use the ${if defined($var) then val1 else val2}
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
            sec_expression = f'{expression} + "{s.replace("^", "")}"'

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
                sec_expression = '${sub({inp}, "\\\\{old_ext}$", "{new_ext}")}'.format(
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
    node2, step_identifier: str, step_alias: str, resource_overrides: Dict[str, str]
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

    ins = node.inputs()

    # Let's sanity check our inputs, make sure they're all present:
    missing_keys = [
        k
        for k in ins.keys()
        if k not in node.sources and not (ins[k].input_type.optional or ins[k].default)
    ]
    if missing_keys:
        raise Exception(
            f"Error when building connections for step '{node.id()}', "
            f"missing the required connection(s): '{', '.join(missing_keys)}'"
        )

    # One step => One WorkflowCall. We need to traverse the edge list to see if there's a scatter
    # then we can build up the WorkflowCall and then wrap them in scatters
    scatterable: List[StepTagInput] = []

    if node.scatter:
        scatterable = [node.sources[k] for k in node.scatter.fields]

        for si in scatterable:
            if si.multiple_inputs or isinstance(si.source(), list):
                raise NotImplementedError(
                    f"The edge '{si.dotted_source()}' on node '{node.id()}' scatters"
                    f"on multiple inputs, and I don't know how this can be implemented in WDL"
                )

    # We need to replace the scatterable key(s) with some random variable, eg: for i in iterable:
    scattered_ordered_variable_identifiers = [
        "i",
        "j",
        "k",
        "x",
        "y",
        "z",
        "a",
        "b",
        "c",
        "ii",
        "jj",
        "kk",
        "xx",
        "yy",
        "zz",
    ]
    scattered_old_to_new_identifier = {
        k.dotted_source(): (k.dotted_source(), k.source())
        for k in scatterable
        if not isinstance(k.dotted_source(), list)
    }
    current_identifiers_that_are_scattered = set(
        v[0] for v in scattered_old_to_new_identifier.values()
    )

    # We'll wrap everything in the scatter block later, but let's replace the fields we need to scatter
    # with the new scatter variable (we'll try to guess one based on the fieldname). We might need to eventually
    # pass the workflow inputs to make sure now conflict will arise.
    # Todo: Pass Workflow input tags to wdl scatter generation to ensure scatter var doesn't conflict with inputs

    for s in scatterable:
        current_identifiers_that_are_scattered.remove(s.dotted_source())
        e: Edge = first_value(s.source_map)
        new_var = e.stag[0] if e.stag else e.source_dotted()[0]

        while new_var in current_identifiers_that_are_scattered:
            new_var = scattered_ordered_variable_identifiers.pop(0)
        scattered_old_to_new_identifier[s.dotted_source()] = (new_var, e.start)
        current_identifiers_that_are_scattered.add(new_var)

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
                    t = (
                        first_value(x.start.outputs())
                        if not x.stag
                        else x.start.outputs()[x.stag]
                    )

                    unique_types.update(t.output_type.secondary_files() or [""])
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
                    f.input_type.subtype().secondary_files()
                    if isinstance(f.input_type, Array)
                    else f.input_type.secondary_files()
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
            it = source.finish.inputs()[source.ftag].input_type
            if source.stag:
                ot = source.start.outputs()[source.stag].output_type
            else:
                ot = first_value(source.start.outputs()).output_type
            if (
                isinstance(it, Array)
                and not isinstance(ot, Array)
                and not source.scatter
            ):
                array_input_from_single_source = True
        secondary = None
        # We're connecting to another step
        if source and isinstance(source.finish, StepNode) and source.ftag:

            it = source.finish.inputs()[source.ftag].input_type

            secondary = (
                it.subtype().secondary_files()
                if isinstance(it, Array)
                else it.secondary_files()
            )
            if secondary and isinstance(source.start, StepNode) and source.stag:

                ot = source.start.outputs()[source.stag].output_type

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
                ] = f"{identifier}[{idx + 1}]"

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

    # So, the current way of mixing accessory files is not really                        _
    # supported, but a little complicated basically, if our scatterable edge           _| |
    # contains secondary files, they'll all be arrays of separate files, eg:         _| | |
    #                                                                               | | | | __
    # File[] bams = [...]                                                           | | | |/  \
    # File[] bais = [...]                                                           |       /\ \
    #                                                                               |       \/ /
    # We can handle this by transposing the array of both items, eg:                 \        /
    #                                                                                 |      /
    #     transpose([bam1, bam2, ..., bamn], [bai1, bai2, ..., bai3])                 |     |
    #           => [[bam1, bai1], [bam2, bai2], ..., [bamn, bain]]
    #
    # and then unwrap them using their indices and hoefully have everything line up:
    #
    # Source: https://software.broadinstitute.org/wdl/documentation/spec#arrayarrayx-transposearrayarrayx

    for s in scatterable:
        if not isinstance(s.finish, StepNode) or not s.ftag:
            raise Exception(
                "An internal error has occured when generating scatterable input map"
            )
        secondary = s.finish.tool.inputs_map()[s.ftag].input_type.secondary_files()
        if secondary:
            ds = s.dotted_source()
            joined_tags = ", ".join(
                get_secondary_tag_from_original_tag(ds, sec) for sec in secondary
            )
            transformed = f"transform([{ds}, {joined_tags}])"
            call = wdl.WorkflowScatter(
                scattered_old_to_new_identifier[s.dotted_source()][0],
                transformed,
                [call],
            )

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

                insource = f"select_first([{insource}, {resolved}])"

            call = wdl.WorkflowScatter(newid, insource, [call])

    return call


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
        return (
            value.generated_filename()
            if string_environment
            else f'"{value.generated_filename()}"'
        )
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

    raise Exception("Could not detect type %s to convert to input value" % type(value))


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
    name = resolve_tool_input_value(inp, **debugkwargs)
    if string_environment:
        return f"${{{name}}}"
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
    tool_name_to_tool: Dict[str, Tool] = {t.id().lower(): t for t in tools}
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


def apply_secondary_file_format_to_filename(
    filepath: Optional[str], secondary_file: str
):
    """
    This is actually clever, you can probably trust this to do what you want.
    :param filepath: Filename to base
    :param secondary_file: CWL secondary format (Remove 1 extension for each leading ^.
    """
    if not filepath:
        return None

    fixed_sec = secondary_file.lstrip("^")
    leading = len(secondary_file) - len(fixed_sec)
    if leading <= 0:
        return filepath + fixed_sec

    basepath = ""
    filename = filepath
    if "/" in filename:
        idx = len(filepath) - filepath[::-1].index("/")
        basepath = filepath[:idx]
        filename = filepath[idx:]

    split = filename.split(".")
    return basepath + ".".join(split[: -min(leading, len(split) - 1)]) + fixed_sec


def prepare_move_statement_for_potential_secondaries(ti: ToolInput):
    it = ti.input_type
    if not issubclass(type(it), File) or not it.secondary_files():
        return []

    commands = []
    for s in it.secondary_files():
        sectag = get_secondary_tag_from_original_tag(ti.id(), s)
        commands.append(
            wdl.Task.Command(
                f'if [ $(dirname "${{{sectag}}}") != $(dirname "{ti.id()}") ]; then mv ${{{sectag}}} $(dirname ${{{ti.id()}}}); fi'
            )
        )
    return commands


def prepare_move_statement_for_input_to_localise(ti: ToolInput):
    it = ti.input_type

    if issubclass(type(it), File):
        commands = [wdl.Task.Command(f"mv ${{{ti.id()}}} -t .")]
        if it.secondary_files():
            for s in it.secondary_files():
                commands.append(
                    wdl.Task.Command(
                        f"mv ${{{get_secondary_tag_from_original_tag(ti.id(), s)}}} -t ."
                    )
                )
        return commands
    if isinstance(it, Array) and issubclass(type(it.subtype()), File):
        subtype = it.subtype()
        commands = [wdl.Task.Command("mv ${{sep=' ' {s}}} -t .".format(s=ti.id()))]
        if subtype.secondary_files():
            for s in it.secondary_files():
                commands.append(
                    wdl.Task.Command(
                        "mv ${{sep=' ' {s}}} -t .".format(
                            s=get_secondary_tag_from_original_tag(ti.id(), s)
                        )
                    )
                )

        return commands

    raise Exception(f"WDL is unable to localise type '{type(it)}'")


def build_resource_override_maps_for_workflow(wf, prefix=None) -> List[wdl.Input]:
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
                    wdl.Input(
                        wdl.WdlType.parse_type("Int?"), tool_pre + "runtime_memory"
                    ),
                    wdl.Input(wdl.WdlType.parse_type("Int?"), tool_pre + "runtime_cpu"),
                    wdl.Input(
                        wdl.WdlType.parse_type("String"), tool_pre + "runtime_disks"
                    ),
                ]
            )
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            inputs.extend(build_resource_override_maps_for_workflow(tool, tool_pre))

    return inputs
