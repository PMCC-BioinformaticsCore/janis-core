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

from janis_core.deps import wdlgen as wdl

from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.operators.logical import If, IsDefined
from janis_core.operators.standard import FirstOperator
from janis_core.types import get_instantiated_type, DataType

from janis_core.types.data_types import is_python_primitive

from janis_core.code.codetool import CodeTool
from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.tool.commandtool import CommandTool, ToolInput, ToolArgument, ToolOutput
from janis_core.tool.tool import Tool, TOutput, ToolType
from janis_core.translations.translationbase import (
    TranslatorBase,
    TranslatorMeta,
    try_catch_translate,
)
from janis_core.operators import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    MemorySelector,
    StringFormatter,
    Selector,
    Operator,
    StepOutputSelector,
    InputNodeSelector,
    TimeSelector,
    DiskSelector,
    ResourceSelector,
    AliasSelector,
)
from janis_core.types.common_data_types import (
    Stdout,
    Stderr,
    Array,
    Boolean,
    Filename,
    File,
    Directory,
    Int,
    Float,
    Double,
    String,
)
from janis_core.utils import (
    first_value,
    recursive_2param_wrap,
    find_duplicates,
    generate_cat_command_from_statements,
)
from janis_core.utils.generators import generate_new_id_from
from janis_core.utils.logger import Logger
from janis_core.utils.pickvalue import PickValue
from janis_core.utils.scatter import ScatterDescription, ScatterMethod
from janis_core.utils.validators import Validators
from janis_core.utils.secondary import (
    split_secondary_file_carats,
    apply_secondary_file_format_to_filename,
)

# from janis_core.tool.step import StepNode


## PRIMARY TRANSLATION METHODS
from janis_core.workflow.workflow import InputNode, StepNode

SED_REMOVE_EXTENSION = "| sed 's/\\.[^.]*$//'"
REMOVE_EXTENSION = (
    lambda x, iterations: f"$(echo '{x}' {iterations * SED_REMOVE_EXTENSION})"
)


class CustomGlob(Selector):
    def __init__(self, expression):
        self.expression = expression

    def returntype(self):
        return Array(File)

    def to_string_formatter(self):
        raise Exception("Not supported for CustomGlob")


class WdlTranslator(TranslatorBase, metaclass=TranslatorMeta):
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
    @try_catch_translate(type="workflow")
    def translate_workflow(
        cls,
        wfi,
        with_container=True,
        with_resource_overrides=False,
        is_nested_tool=False,
        allow_empty_container=False,
        container_override=None,
    ) -> Tuple[wdl.Workflow, Dict[str, any]]:
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

        inputsdict = {t.id(): ToolInput(t.id(), t.intype) for t in wf.tool_inputs()}

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

            if i.default is not None:
                expr = WdlTranslator.unwrap_expression(i.default, inputsdict=inputsdict)

            wd = dt.wdl(has_default=i.default is not None)

            w.inputs.append(
                wdl.Input(
                    data_type=wd, name=i.id(), expression=expr, requires_quotes=False
                )
            )

            is_array = i.datatype.is_array()
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
            outtag = None

            if isinstance(o.source, list):
                outtag = WdlTranslator.unwrap_expression(
                    expression=o.source,
                    inputsdict=None,
                    string_environment=False,
                    wfid=wf.id(),
                    outputid=o.id(),
                )

            else:
                outtag = WdlTranslator.unwrap_expression(
                    expression=o.source,
                    inputsdict=None,
                    string_environment=False,
                    wfid=wf.id(),
                    outputid=o.id(),
                )

            w.outputs.append(wdl.Output(o.datatype.wdl(), o.id(), outtag))

            fundamental_outtype = o.datatype
            if fundamental_outtype.is_array():
                fundamental_outtype = fundamental_outtype.fundamental_type()
            if fundamental_outtype.secondary_files():
                if isinstance(o.source, InputNodeSelector):
                    src = [o.source.id()]
                elif isinstance(o.source, StepOutputSelector):
                    src = [o.source.node.id(), o.source.tag]
                else:
                    raise Exception(
                        f"Unsupported type for output with secondary files: {type(o.source)}"
                    )
                w.outputs.extend(
                    wdl.Output(
                        o.datatype.wdl(),
                        get_secondary_tag_from_original_tag(o.id(), s),
                        ".".join(
                            [*src[:-1], get_secondary_tag_from_original_tag(src[-1], s)]
                        ),
                    )
                    for s in fundamental_outtype.secondary_files()
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
                if t.type() == ToolType.Workflow:
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
                node2=s,
                step_identifier=tool_aliases[t.versioned_id().lower()].upper()
                + "."
                + t.id(),
                resource_overrides=resource_overrides,
                invalid_identifiers=forbiddenidentifiers,
                inputsdict=inputsdict,
            )

            w.calls.append(call)

        return w, wtools

    @classmethod
    @try_catch_translate(type="command tool")
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
        outs: List[wdl.Output] = cls.translate_tool_outputs(toolouts, inmap, tool=tool)
        command_args = cls.translate_tool_args(
            tool.arguments(), inmap, tool=tool, toolId=tool.id()
        )
        command_ins = cls.build_command_from_inputs(tool.inputs())

        commands = [wdl.Task.Command("set -e")]

        # generate directories / file to create commands
        commands.extend(cls.build_commands_for_file_to_create(tool))

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

        namedwdlouts = {t.name: t for t in outs}
        for to in toolouts:
            commands.extend(
                prepare_move_statements_for_output(to, namedwdlouts[to.id()].expression)
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
        cls.add_runtimefield_overrides_for_wdl(
            runtime_block=r,
            tool=tool,
            inmap=inmap,
            with_resource_overrides=with_resource_overrides,
        )

        return wdl.Task(tool.id(), ins, outs, commands, r, version="development")

    @classmethod
    @try_catch_translate(type="code tool")
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

        ins = cls.get_resource_override_inputs() + [
            ToolInput(
                t.id(),
                input_type=t.intype,
                prefix=f"--{t.id()}",
                default=t.default,
                doc=t.doc,
            )
            for t in tool.tool_inputs()
        ]

        tr_ins = cls.translate_tool_inputs(ins)

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

        tr_outs = cls.translate_tool_outputs(outs, {}, tool=tool)

        commands = []

        scriptname = tool.script_name()

        commands.append(
            wdl.Task.Command(
                f"""
cat <<EOT >> {scriptname}
    {tool.prepared_script(SupportedTranslation.WDL)}
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

        inmap = {t.id(): t for t in ins}
        cls.add_runtimefield_overrides_for_wdl(
            r, tool=tool, inmap=inmap, with_resource_overrides=with_resource_overrides
        )

        return wdl.Task(tool.id(), tr_ins, tr_outs, commands, r, version="development")

    @staticmethod
    def wrap_if_string_environment(value, string_environment: bool):
        return f'"{value}"' if not string_environment else value

    @classmethod
    def unwrap_expression(
        cls,
        expression,
        inputsdict=None,
        string_environment=False,
        tool=None,
        for_output=False,
        **debugkwargs,
    ):
        if expression is None:
            return ""

        wrap_in_code_block = lambda x: f"~{{{x}}}" if string_environment else x

        if isinstance(expression, StepNode):
            raise Exception(
                f"The Step node '{expression.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )

        if isinstance(expression, list):
            toolid = value_or_default(debugkwargs.get("tool_id"), "get-value-list")
            joined_values = ", ".join(
                str(
                    cls.unwrap_expression(
                        expression[i],
                        inputsdict,
                        string_environment=False,
                        tool_id=toolid + "." + str(i),
                    )
                )
                for i in range(len(expression))
            )
            return f"[{joined_values}]"
        if is_python_primitive(expression):
            if isinstance(expression, str):
                if string_environment:
                    return expression
                return cls.wrap_if_string_environment(
                    prepare_escaped_string(expression), string_environment
                )
            if isinstance(expression, bool):
                return "true" if expression else "false"

            return str(expression)
        elif isinstance(expression, Filename):
            gen_filename = expression.generated_filename(
                replacements={
                    "prefix": WdlTranslator.unwrap_expression(
                        expression.prefix,
                        inputsdict=inputsdict,
                        string_environment=True,
                        for_output=for_output,
                    )
                }
            )
            return cls.wrap_if_string_environment(gen_filename, string_environment)

        elif isinstance(expression, AliasSelector):
            return cls.unwrap_expression(
                expression.inner_selector,
                string_environment=string_environment,
                inputsdict=inputsdict,
                tool=tool,
                for_output=for_output,
                **debugkwargs,
            )

        elif isinstance(expression, StringFormatter):
            return translate_string_formatter(
                selector=expression,
                inputsdict=inputsdict,
                string_environment=string_environment,
                tool=tool,
                **debugkwargs,
            )
        elif isinstance(expression, WildcardSelector):
            raise Exception(
                f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
            )

        elif isinstance(expression, ResourceSelector):

            if not tool:
                raise Exception(
                    f"Tool must be provided when unwrapping ResourceSelector: {type(expression).__name__}"
                )
            operation = expression.get_operation(tool, hints={})
            return cls.unwrap_expression(
                operation,
                string_environment=string_environment,
                inputsdict=inputsdict,
                tool=tool,
                for_output=for_output,
                **debugkwargs,
            )

        elif isinstance(expression, InputSelector):
            if for_output:
                val = prepare_filename_replacements_for(
                    expression, inputsdict=inputsdict
                )
                return wrap_in_code_block(val)
            return translate_input_selector(
                selector=expression,
                inputsdict=inputsdict,
                string_environment=string_environment,
                **debugkwargs,
            )
        elif callable(getattr(expression, "wdl", None)):
            return expression.wdl()

        unwrap_expression_wrap = lambda exp: cls.unwrap_expression(
            exp,
            inputsdict,
            string_environment=False,
            tool=tool,
            for_output=for_output,
            **debugkwargs,
        )

        if isinstance(expression, InputNodeSelector):
            value = expression.input_node.id()
            if expression.input_node.default is not None:
                unwrapped_default = unwrap_expression_wrap(
                    expression.input_node.default
                )
                value = f"select_first([{value}, {unwrapped_default}])"
            return wrap_in_code_block(value)

        if isinstance(expression, StepOutputSelector):
            value = expression.node.id() + "." + expression.tag
            return wrap_in_code_block(value)

        elif isinstance(expression, Operator):
            return wrap_in_code_block(
                expression.to_wdl(unwrap_expression_wrap, *expression.args)
            )

        warning = ""
        if isclass(expression):
            stype = expression.__name__
            warning = f", this is likely due to the '{stype}' not being initialised"
        else:
            stype = expression.__class__.__name__
        raise Exception(
            f"Could not convert expression '{expression}' as could detect type '{stype}' to convert to input value{warning}"
        )

    @classmethod
    def unwrap_expression_for_output(
        cls,
        output: ToolOutput,
        expression,
        inputsdict=None,
        string_environment=False,
        **debugkwargs,
    ):
        """
        :param output:
        :param expression:
        :param inputsdict:
        :param string_environment:
        :param debugkwargs:
        :return:
        """
        if isinstance(expression, CustomGlob):
            return expression.expression
        elif isinstance(output.output_type, Stdout) or isinstance(expression, Stdout):
            # can't get here with secondary_format
            return "stdout()"
        elif isinstance(output.output_type, Stderr) or isinstance(expression, Stderr):
            return "stderr()"

        if isinstance(expression, list):
            toolid = value_or_default(debugkwargs.get("tool_id"), "get-value-list")
            joined_values = ", ".join(
                str(
                    cls.unwrap_expression_for_output(
                        output=output,
                        expression=expression[i],
                        inputsdict=inputsdict,
                        string_environment=False,
                        tool_id=toolid + "." + str(i),
                    )
                )
                for i in range(len(expression))
            )
            return f"[{joined_values}]"
        if is_python_primitive(expression):
            if isinstance(expression, str):
                return cls.wrap_if_string_environment(expression, string_environment)
            if isinstance(expression, bool):
                return "true" if expression else "false"

            return str(expression)
        elif isinstance(expression, StringFormatter):
            return translate_string_formatter_for_output(
                out=output,
                selector=expression,
                inputsdict=inputsdict,
                string_environment=string_environment,
                **debugkwargs,
            )
        elif isinstance(expression, WildcardSelector):
            is_single = not output.output_type.is_array()
            select_first = None
            is_single_optional = None
            if is_single:
                is_single_optional = output.output_type.optional
                if not expression.select_first:
                    Logger.info(
                        f"The command tool ({debugkwargs}).{output.tag}' used a star-bind (*) glob to find the output, "
                        f"but the return type was not an array. For WDL, the first element will be used, "
                        f"ie: 'glob(\"{expression.wildcard}\")[0]'"
                    )
                    select_first = True

            base_expression = translate_wildcard_selector(
                expression,
                override_select_first=select_first,
                is_optional=is_single_optional,
            )

            return base_expression

        elif isinstance(expression, InputSelector):
            return translate_input_selector_for_output(
                out=output,
                selector=expression,
                inputsdict=inputsdict,
                string_environment=string_environment,
                **debugkwargs,
            )
        elif callable(getattr(expression, "wdl", None)):
            return expression.wdl()

        wrap_in_code_block = lambda x: f"~{{{x}}}" if string_environment else x
        unwrap_expression_wrap = lambda exp: cls.unwrap_expression_for_output(
            output=output,
            expression=exp,
            inputsdict=inputsdict,
            string_environment=False,
            **debugkwargs,
        )

        if isinstance(expression, (StepOutputSelector, InputNodeSelector)):
            raise Exception(
                "An InputnodeSelector or StepOutputSelector cannot be used to glob outputs"
            )

        elif isinstance(expression, Operator):
            return wrap_in_code_block(
                expression.to_wdl(unwrap_expression_wrap, *expression.args)
            )

        warning = ""
        if isclass(expression):
            stype = expression.__name__
            warning = f", this is likely due to the '{stype}' not being initialised"
        else:
            stype = expression.__class__.__name__

        raise Exception(
            f"Tool ({debugkwargs}) has an unrecognised glob type: '{stype}' ({expression}), this is "
            f"deprecated. Please use the a Selector to build the outputs for '{output.id()}'"
            + warning
        )

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
                    if i.input_type.is_array()
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
        cls, tooloutputs: List[ToolOutput], inputsmap: Dict[str, ToolInput], tool
    ):
        outs: List[wdl.Output] = []

        for o in tooloutputs:
            wdl_type = wdl.WdlType.parse_type(o.output_type.wdl())
            expression = cls.unwrap_expression_for_output(
                o, o.selector, inputsdict=inputsmap, tool=tool, toolid=tool.id()
            )
            outs.append(wdl.Output(wdl_type, o.id(), expression))
            outs.extend(
                cls.prepare_secondary_tool_outputs(
                    out=o,
                    original_expression=o.selector,
                    expression=expression,
                    toolid=tool.id(),
                )
            )

        return outs

    @classmethod
    def prepare_secondary_tool_outputs(
        cls, out: ToolOutput, original_expression: any, expression: str, toolid: str
    ) -> List[wdl.Output]:
        if not (
            isinstance(out.output_type, File) and out.output_type.secondary_files()
        ):
            return []

        islist = isinstance(expression, list)

        if (
            out.output_type.is_array()
            and isinstance(out.output_type.subtype(), File)
            and out.output_type.subtype().secondary_files()
        ):
            if isinstance(original_expression, WildcardSelector):
                # do custom override for wildcard selector
                is_single = not out.output_type.is_array()
                select_first = None
                is_single_optional = None
                if is_single and not original_expression.select_first:
                    select_first = True
                    is_single_optional = out.output_type.optional
                ftype = out.output_type.subtype().wdl()
                return [
                    wdl.Output(
                        ftype,
                        get_secondary_tag_from_original_tag(out.id(), s),
                        translate_wildcard_selector(
                            original_expression,
                            secondary_format=s,
                            override_select_first=select_first,
                            is_optional=is_single_optional,
                        ),
                    )
                    for s in out.output_type.subtype().secondary_files()
                ]
            elif islist:
                Logger.info(
                    "Special handling for an Array return type with a list expressions"
                )
            else:
                raise Exception(
                    "Janis isn't sure how to collect secondary files for an array yet"
                )

        outs = []
        if isinstance(out.output_type, File) and out.output_type.secondary_files():
            # eep we have secondary files
            ot = get_instantiated_type(out.output_type)
            ftype = ot.wdl()
            for s in ot.secondary_files():
                tag = get_secondary_tag_from_original_tag(out.id(), s)
                ar_exp = expression if islist else [expression]
                potential_extensions = ot.get_extensions()
                if "^" not in s:
                    exp = [(ex + f' + "{s}"') for ex in ar_exp]
                elif potential_extensions:
                    exp = []
                    for ex in ar_exp:
                        inner_exp = ex
                        for ext in potential_extensions:
                            inner_exp = 'sub({inp}, "\\\\{old_ext}$", "{new_ext}")'.format(
                                inp=inner_exp, old_ext=ext, new_ext=s.replace("^", "")
                            )
                        exp.append(inner_exp)

                else:
                    raise Exception(
                        f"Unsure how to handle secondary file '{s}' for the tool output '{out.id()}' (ToolId={toolid})"
                        f" as it uses the escape characater '^' but Janis can't determine the extension of the output."
                        f"This could be resolved by ensuring the definition for '{ot.__class__.__name__}' contains an extension."
                    )

                if ot.optional:
                    exp = [
                        f"if defined({expression}) then ({e}) else None" for e in exp
                    ]

                outs.append(wdl.Output(ftype, tag, exp if islist else exp[0]))

        return outs

    @classmethod
    def translate_tool_args(
        cls,
        toolargs: List[ToolArgument],
        inpmap: Dict[str, ToolInput],
        tool,
        **debugkwargs,
    ):
        if not toolargs:
            return []
        commandargs = []
        for a in toolargs:
            val = cls.unwrap_expression(
                a.value, inpmap, tool=tool, string_environment=True
            )
            should_wrap_in_quotes = isinstance(val, str) and (
                a.shell_quote is None or a.shell_quote
            )
            wrapped_val = f"'{val}'" if should_wrap_in_quotes else val
            commandargs.append(
                wdl.Task.Command.CommandArgument.from_fields(
                    a.prefix, wrapped_val, a.position
                )
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
    def build_commands_for_file_to_create(
        cls, tool: CommandTool
    ) -> List[wdl.Task.Command]:
        commands = []
        inputsdict = {t.id(): t for t in tool.inputs()}

        directories = tool.directories_to_create()
        files = tool.files_to_create()

        if directories is not None:
            directories = (
                directories if isinstance(directories, list) else [directories]
            )
            for directory in directories:
                unwrapped_dir = cls.unwrap_expression(
                    directory, inputsdict=inputsdict, tool=tool, string_environment=True
                )
                commands.append(f"mkdir -p '{unwrapped_dir}'")
        if files:
            for path, contents in files if isinstance(files, list) else files.items():
                unwrapped_path = cls.unwrap_expression(
                    path, inputsdict=inputsdict, tool=tool, string_environment=True
                )
                unwrapped_contents = cls.unwrap_expression(
                    contents, inputsdict=inputsdict, tool=tool, string_environment=True
                )
                commands.append(
                    generate_cat_command_from_statements(
                        path=unwrapped_path, contents=unwrapped_contents
                    )
                )

        return list(map(wdl.Task.Command, commands))

    @classmethod
    def add_runtimefield_overrides_for_wdl(
        cls, runtime_block, tool, inmap, with_resource_overrides
    ):
        runtime_block.kwargs["cpu"] = cls.unwrap_expression(
            CpuSelector(), inmap, string_environment=False, tool=tool, id="runtimestats"
        )
        runtime_block.kwargs["memory"] = cls.unwrap_expression(
            StringFormatter("{value}G", value=MemorySelector()),
            inmap,
            string_environment=False,
            tool=tool,
            id="runtimestats",
        )
        runtime_block.kwargs["duration"] = cls.unwrap_expression(
            TimeSelector(),
            inmap,
            string_environment=False,
            tool=tool,
            id="runtimestats",
        )
        runtime_block.kwargs["disks"] = cls.unwrap_expression(
            StringFormatter("local-disk {value} SSD", value=DiskSelector()),
            inmap,
            string_environment=False,
            tool=tool,
            id="runtimestats",
        )

        if with_resource_overrides:
            runtime_block.kwargs["zones"] = '"australia-southeast1-b"'

        runtime_block.kwargs["preemptible"] = 2

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
        max_duration=None,
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
        is_workflow = tool.type() == ToolType.Workflow

        if is_workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value is not None
                or (i.default is not None and not isinstance(i.default, Selector))
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
            elif i.intype.is_array() and i.intype.subtype().secondary_files():
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
                    tool,
                    hints,
                    max_cores=max_cores,
                    max_mem=max_mem,
                    max_duration=max_duration,
                    inputs=ad,
                    is_root=True,
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
        max_duration=None,
        inputs=None,
        prefix=None,
        is_root=False,
    ):
        from janis_core.workflow.workflow import Workflow

        is_workflow = tool.type() == ToolType.Workflow
        d = super().build_resources_input(
            tool=tool,
            hints=hints,
            max_cores=max_cores,
            max_mem=max_mem,
            max_duration=max_duration,
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


def prepare_escaped_string(value: str):
    return json.dumps(value)[1:-1]


def resolve_tool_input_value(
    tool_input: ToolInput, inputsdict, string_environment=False, **debugkwargs
):
    name = tool_input.id()
    indefault = (
        tool_input.input_type
        if isinstance(tool_input.input_type, Filename)
        else tool_input.default
    )

    default = None
    if isinstance(indefault, ResourceSelector):

        if indefault.default:
            default = (
                f"select_first([{indefault.input_to_select}, {str(indefault.default)}])"
            )
        else:
            default = indefault.input_to_select

    elif isinstance(indefault, InputSelector):
        Logger.critical(
            f"WDL does not support command line level defaults that select a different input, this will remove the "
            f"value: '{indefault}' for tool_input '{tool_input.tag}'"
        )

    elif indefault is not None:
        default = WdlTranslator.unwrap_expression(
            indefault,
            inputsdict=inputsdict,
            string_environment=string_environment,
            **debugkwargs,
        )

    if default is not None:
        # Default should imply optional input
        name = f"select_first([{name}, {default}])"

    if tool_input.localise_file:
        if tool_input.input_type.is_array():
            raise Exception(
                "Localising files through `basename(x)` is unavailable for arrays of files: https://github.com/openwdl/wdl/issues/333"
            )
        if tool_input.presents_as:
            return (
                tool_input.presents_as
                if string_environment
                else f'"{tool_input.presents_as}"'
            )
        else:
            name = "basename(%s)" % name

    return f"~{{{name}}}" if string_environment else name


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_input.position is not None or tool_input.prefix):
        return None

    if tool_input.localise_file and tool_input.presents_as:
        return wdl.Task.Command.CommandInput(
            value=tool_input.presents_as, position=tool_input.position
        )

    name = resolve_tool_input_value(
        tool_input, inputsdict=inputsdict, string_environment=False, **debugkwargs
    )
    intype = tool_input.input_type

    optional = (not isinstance(intype, Filename) and intype.optional) or (
        isinstance(tool_input.default, CpuSelector) and tool_input.default is None
    )
    position = tool_input.position

    separate_value_from_prefix = tool_input.separate_value_from_prefix is not False
    prefix = tool_input.prefix if tool_input.prefix else ""
    tprefix = prefix

    intype = tool_input.input_type

    is_flag = isinstance(intype, Boolean)

    if prefix and separate_value_from_prefix and not is_flag:
        tprefix += " "

    expr = name

    if isinstance(intype, Boolean):
        if tool_input.prefix:
            if tool_input.default is not None or not intype.optional:
                condition = expr
            else:
                condition = f"(defined({expr}) && select_first([{expr}]))"
            expr = f'~{{if {condition} then "{tprefix}" else ""}}'
    elif intype.is_array():

        separator = tool_input.separator if tool_input.separator is not None else " "
        should_quote = (
            isinstance(intype.subtype(), (String, File, Directory))
            and tool_input.shell_quote is not False
        )
        condition_for_binding = None

        if intype.optional:
            rexpr = expr
            expr = f"select_first([{expr}])"
            condition_for_binding = f"(defined({rexpr}) && length({expr}) > 0)"
        else:
            condition_for_binding = f"length({expr}) > 0"

        if intype.subtype().optional:
            expr = f"select_all({expr})"

        if should_quote:
            if tool_input.prefix_applies_to_all_elements:
                separator = f"'{separator}{tprefix}'"
            else:
                separator = f"'{separator}'"

            if tprefix:
                expr = f'"{tprefix}\'" + sep("{separator}", {expr}) + "\'"'
            else:
                expr = f'"\'" + sep("{separator}", {expr}) + "\'"'

        else:
            if tprefix:
                if tool_input.prefix_applies_to_all_elements:
                    expr = f'sep("{separator}", prefix("{tprefix}", {expr}))'
                else:
                    expr = f'"{tprefix}" + sep("{separator}", {expr})'
            else:
                expr = f'sep("{separator}", {expr})'

        if condition_for_binding is not None:
            expr = f'~{{if {condition_for_binding} then {expr} else ""}}'
        else:
            expr = f"~{{{expr}}}"
    elif (
        isinstance(intype, (String, File, Directory))
        and tool_input.shell_quote is not False
    ):
        if tprefix:
            if optional:
                expr = f'~{{if defined({expr}) then ("{tprefix}\'" + {expr} + "\'") else ""}}'
            else:
                expr = f"{tprefix}'~{{{expr}}}'"
        else:
            if optional:
                expr = f'~{{if defined({expr}) then ("\'" + {expr} + "\'") else ""}}'
            else:
                expr = f"'~{{{expr}}}'"

    else:
        if prefix:
            if optional:
                expr = f"~{{if defined({expr}) then (\"{tprefix}\" + {expr}) else ''}}"
            else:
                expr = f"{tprefix}~{{{expr}}}"
        else:
            expr = f"~{{{expr}}}"

    # there used to be a whole lot of login in the wdl.Task.Command.CommandInput but it's been moved to here now
    return wdl.Task.Command.CommandInput(value=expr, position=tool_input.position)


def translate_input_selector_for_output(
    out: ToolOutput,
    selector: InputSelector,
    inputsdict: Dict[str, ToolInput],
    string_environment=False,
    **debugkwargs,
) -> List[wdl.Output]:
    expression = translate_input_selector(
        selector, inputsdict, string_environment=False, **debugkwargs
    )

    tool_in = inputsdict.get(selector.input_to_select)
    if not tool_in:
        raise Exception(
            f"The InputSelector for tool '{debugkwargs}.{out.id()}' did not select an input (tried: '{selector.input_to_select}')"
        )

    return expression


def translate_input_selector_for_secondary_output(
    out: ToolOutput,
    selector: InputSelector,
    inputsdict: Dict[str, ToolInput],
    **debugkwargs,
) -> List[wdl.Output]:
    expression = translate_input_selector(
        selector, inputsdict, string_environment=False, **debugkwargs
    )

    tool_in = inputsdict.get(selector.input_to_select)
    if not tool_in:
        raise Exception(
            f"The InputSelector for tool '{debugkwargs}.{out.id()}' did not select an input (tried: '{selector.input_to_select}')"
        )

    return expression


def translate_string_formatter_for_output(
    out,
    selector: StringFormatter,
    inputsdict: Dict[str, ToolInput],
    string_environment,
    **debugkwargs,
) -> str:
    """
    The output glob was a string formatter, so we'll need to build the correct glob
    by resolving the string formatter. Some minor complications involve how an output
    with a secondary file must resolve input selectors.

    For example, if you are generating an output with the Filename class, and your output
    type has secondary files, you will need to translate the generated filename with
    respect to the secondary file extension. Or the File class also has a recommended
    "extension" property now that this should consider.

    :param inputsdict:
    :param out:
    :param selector:
    :return:
    """
    inputs_to_retranslate = {
        k: v
        for k, v in selector.kwargs.items()
        if not any(isinstance(v, t) for t in StringFormatter.resolved_types)
    }

    resolved_kwargs = {
        **selector.kwargs,
        **{
            k: WdlTranslator.unwrap_expression(
                v, inputsdict=inputsdict, string_environment=True, **debugkwargs
            )
            for k, v in inputs_to_retranslate.items()
        },
    }

    return f'"{selector.resolve_with_resolved_values(**resolved_kwargs)}"'


def validate_step_with_multiple_sources(node, edge, k, input_name_maps):
    multiple_sources_failure_reasons = []

    unique_types = set()
    for x in edge.source():
        t: DataType = x.source.returntype()
        unique_types.update(t.secondary_files() or [])

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


def translate_step_node(
    node2,
    step_identifier: str,
    resource_overrides: Dict[str, str],
    invalid_identifiers: Set[str],
    inputsdict: Dict[str, any],
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
        if k not in node.sources
        and not (ins[k].intype.optional is True or ins[k].default is not None)
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
            if si.multiple_inputs
            or (isinstance(si.source(), list) and len(si.source()) > 1)
        ]
        if len(invalid_sources) > 0:
            invalid_sources_str = ", ".join(f"{si.source()}" for si in invalid_sources)
            raise NotImplementedError(
                f"The edge(s) '{invalid_sources_str}' on node '{node.id()}' scatters"
                f"on multiple inputs, this behaviour has not been implemented"
            )

    # 1. Generate replacement of the scatterable key(s) with some random variable, eg: for i in iterable:
    #
    #       - Currently, Janis does not support operating on the object to scatter, and there's no mechanism from
    #           operating on the scattered value. See the following GH comment for more information:
    #           (https://github.com/PMCC-BioinformaticsCore/janis-core/pull/10#issuecomment-605807815)
    #

    scattered_old_to_new_identifier = generate_scatterable_details(
        scatterable, forbiddenidentifiers=invalid_identifiers
    )

    # Let's map the inputs, to the source. We're using a dictionary for the map atm, but WDL requires the _format:
    #       fieldName: sourceCall.Output

    inputs_map = {}
    for k, inp in ins.items():
        if k not in node.sources:
            continue

        steptag_input: StepTagInput = node.sources[k]
        intype = inp.intype
        src: Edge = steptag_input.source()  # potentially single item or array

        ar_source = src if isinstance(src, list) else [src]
        # these two are very closely tied, they'll determine whether our
        # input to the step connection is single or an array
        has_multiple_sources = isinstance(src, list) and len(src) > 1
        array_input_from_single_source = False

        if has_multiple_sources:
            # let's do some checks, make sure we're okay
            validate_step_with_multiple_sources(node, steptag_input, k, inputs_map)

        elif ar_source:
            source = ar_source[0]

            ot = source.source.returntype()
            if intype.is_array() and not ot.is_array() and not source.scatter:
                array_input_from_single_source = True
        else:
            Logger.critical(
                f"Skipping connection to '{steptag_input.finish}.{steptag_input.ftag}' had no source or default, "
                f"please raise an issue as investigation may be required"
            )
            continue

        # Checks over, let's continue!

        secondaries = (
            intype.secondary_files()
            if not intype.is_array()
            else intype.subtype().secondary_files()
        ) or []
        # place to put the processed_sources:
        #   Key=None is for the regular input
        #   Key=$sec_tag is for each secondary file
        unwrapped_sources = {k: [] for k in [None, *secondaries]}

        unwrap_helper = lambda exprsn: WdlTranslator.unwrap_expression(
            exprsn,
            inputsdict=inputsdict,
            string_environment=False,
            stepid=step_identifier,
        )

        for edge in ar_source:
            # we have an expression we need to unwrap,
            # it's going to the step_input [k]

            if secondaries:
                ot = edge.source.returntype()

                sec_out = set(
                    value_or_default(
                        ot.subtype().secondary_files()
                        if ot.is_array()
                        else ot.secondary_files(),
                        default=[],
                    )
                )
                sec_in = set(secondaries)
                if not sec_in.issubset(sec_out):
                    raise Exception(
                        f"An error occurred when connecting '{edge.source}' to "
                        f"'{edge.finish.id()}.{edge.ftag}', there were secondary files in the final node "
                        f"that weren't present in the source: {', '.join(sec_out.difference(sec_in))}"
                    )

            unwrapped_exp = unwrap_helper(edge.source)

            default = None
            if isinstance(edge.source, InputNodeSelector):
                default = unwrap_helper(edge.source.input_node.default)

            is_scattered = unwrapped_exp in scattered_old_to_new_identifier

            if is_scattered:
                unwrapped_exp = scattered_old_to_new_identifier[unwrapped_exp][0]

                for idx in range(len(secondaries)):
                    # we restrict that files with secondaries can't be operated on in the step input
                    sec = secondaries[idx]
                    unwrapped_sources[sec].append(f"{unwrapped_exp}[{idx + 1}]")

                if secondaries:
                    unwrapped_exp += "[0]"
            else:
                for sec in secondaries:
                    unwrapped_sources[sec].append(
                        get_secondary_tag_from_original_tag(unwrapped_exp, sec)
                    )

            unwrapped_sources[None].append(unwrapped_exp)

        should_select_first_element = not (
            array_input_from_single_source or has_multiple_sources
        )
        for tag, value in unwrapped_sources.items():
            if tag is None:
                tag = k
            else:
                tag = get_secondary_tag_from_original_tag(k, tag)

            inputs_map[tag] = (
                value[0]
                if should_select_first_element
                else "[" + ", ".join(value) + "]"
            )

    inputs_map.update(resource_overrides)

    call = wdl.WorkflowCall(step_identifier, step_alias, inputs_map)

    if len(scatterable) > 0:
        call = wrap_scatter_call(
            call, node.scatter, scatterable, scattered_old_to_new_identifier
        )

    if node.when is not None:
        condition = WdlTranslator.unwrap_expression(
            node.when, inputsdict=inputsdict, string_environment=False
        )
        call = wdl.WorkflowConditional(condition, [call])
        # determine how to unwrap when

    return call


def generate_scatterable_details(
    scatterable: List[StepTagInput], forbiddenidentifiers: Set[str]
):
    if not scatterable:
        return {}

    # get the reference from a InputNodeSelector or StepOutputSelector
    get_source = lambda e: WdlTranslator.unwrap_expression(e.source)

    # this dictionary is what we're going to use to map our current
    # identifier to the scattered identifier. This step is just the
    # setup, and in the next for loop, we'll
    scattered_old_to_new_identifier = {}
    for k in scatterable:
        srcs = k.source()
        for edge in srcs if isinstance(srcs, list) else [srcs]:
            src = get_source(edge)
            scattered_old_to_new_identifier[src] = (src, edge.source)

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
            scattered_old_to_new_identifier[get_source(s.source())] = (
                newid,
                s.source_map[0].source,
            )
            forbiddenidentifierscopy.add(newid)
    else:

        for s in scatterable:

            # We asserted earlier that the source_map only has one value (through multipleInputs)
            e: Edge = s.source_map[0]

            # if isinstance(e.source, Operator):
            #     raise Exception(
            #         "Currently, Janis doesn't support operating on a value to be scattered"
            #     )

            original_expr = WdlTranslator.unwrap_expression(s.source().source)
            newid = generate_new_id_from(original_expr, forbiddenidentifierscopy)
            evaluated_operator = WdlTranslator.unwrap_expression(
                e.source, string_environment=False
            )
            scattered_old_to_new_identifier[evaluated_operator] = (newid, e)
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
    get_source = lambda e: WdlTranslator.unwrap_expression(e.source)

    insource_ar = []
    for s in scatterable:
        secondary = s.finish.tool.inputs_map()[s.ftag].intype.secondary_files()
        if secondary:
            ds = get_source(s.source())
            joined_tags = ", ".join(
                get_secondary_tag_from_original_tag(ds, sec) for sec in secondary
            )
            transformed = f"transpose([{ds}, {joined_tags}])"
            insource_ar.append(transformed)

        else:
            (newid, startnode) = scattered_old_to_new_identifier[get_source(s.source())]
            insource = get_source(s.source())
            if isinstance(startnode, InputNode) and startnode.default is not None:
                resolved = WdlTranslator.unwrap_expression(
                    startnode.default, scatterstep=insource
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
        method = "zip" if scatter.method == ScatterMethod.dot else "cross"
        insource = recursive_2param_wrap(method, insource_ar)
        alias = scattered_old_to_new_identifier["-"]

    return wdl.WorkflowScatter(alias, insource, [call])


## SELECTOR HELPERS


def translate_string_formatter(
    selector: StringFormatter, inputsdict, string_environment, **debugkwargs
):
    # we should raise an Exception if any of our inputs are optional without a default

    invalid_select_inputs = [
        (k, selector.kwargs[k].input_to_select)
        for k in selector.kwargs
        # Our selector is getting an input
        if isinstance(selector.kwargs[k], InputSelector)
        and not isinstance(
            selector.kwargs[k],
            (CpuSelector, MemorySelector, DiskSelector, TimeSelector),
        )
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
        Logger.log(
            f'There might be an error when resolving the format "{selector._format}", the tag(s) {tags} respectively '
            f"selected input(s) {inps} that were optional and did NOT have a default value. This might be okay if "
            f"{tags} was wrapped in a IfDefined operator"
        )

    value = selector.resolve_with_resolved_values(
        **{
            k: WdlTranslator.unwrap_expression(
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

    intype = inp.input_type
    if selector.remove_file_extension and (
        File().can_receive_from(intype) or Directory().can_receive_from(intype)
    ):
        if isinstance(intype, File):
            extensions = {
                e
                for e in [intype.extension, *(intype.alternate_extensions or [])]
                if e is not None
            }
            if extensions:
                for ext in extensions:
                    name = f'basename({name}, "{ext}")'
            else:
                name = f"basename({name})"
        else:
            name = f"basename({name})"

    if string_environment:

        return f"~{{{name}}}"
    else:
        return name


# These are now handled entirely by the translate_input_selector_class

# def translate_cpu_selector(selector: CpuSelector, inputsdict, string_environment=True):
#     return translate_input_selector(selector, inputsdict, string_environment=string_environment)
# def translate_mem_selector(selector: MemorySelector, inputsdict, string_environment=True):
#     return translate_input_selector(selector, inputsdict, string_environment=string_environment)


def translate_wildcard_selector(
    selector: WildcardSelector,
    secondary_format: Optional[str] = None,
    override_select_first: Optional[bool] = None,
    is_optional: Optional[bool] = None,
):
    if not selector.wildcard:
        raise Exception(
            "No wildcard was provided for wildcard selector: " + str(selector)
        )

    wildcard = selector.wildcard
    if secondary_format:
        wildcard = apply_secondary_file_format_to_filename(wildcard, secondary_format)

    gl = f'glob("{wildcard}")'
    if selector.select_first or override_select_first:
        if is_optional:
            gl = f"if length({gl}) > 0 then {gl}[0] else None"
        else:
            gl += "[0]"

    return gl


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


def get_secondary_tag_from_original_tag(original, secondary) -> str:
    secondary_without_punctuation = secondary.replace(".", "").replace("^", "")
    return original + "_" + secondary_without_punctuation


def prepare_env_var_setters(
    reqs: Dict[str, Any], inputsdict, **debugkwargs
) -> List[wdl.Task.Command]:
    if not reqs:
        return []

    statements = []
    for k, v in reqs.items():
        val = WdlTranslator.unwrap_expression(
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
            newlocation = f"`dirname ~{{{ti.id()}}}`/{ti.presents_as}"
            base = newlocation
        else:
            newlocation = ti.presents_as
            base = ti.presents_as

        commands.append(wdl.Task.Command(f"cp -f '~{{{ti.id()}}}' '{newlocation}'"))

    if it.secondary_files():
        sec_presents_as = ti.secondaries_present_as or {}

        for s in it.secondary_files():
            sectag = get_secondary_tag_from_original_tag(ti.id(), s)
            if ti.localise_file and not ti.presents_as:
                # move into the current directory
                dest = "."
            else:
                newext, iters = split_secondary_file_carats(sec_presents_as.get(s, s))
                dest = REMOVE_EXTENSION(base, iters) + newext

            commands.append(wdl.Task.Command(f"cp -f '~{{{sectag}}}' {dest}"))

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

            commands.append(
                wdl.Task.Command(
                    f"if [ -f {oldpath} ]; then ln -f {oldpath} {newpath}; fi"
                )
            )

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
                wdl.Input(wdl.WdlType.parse_type("Int?"), prefix + "runtime_disks"),
                wdl.Input(wdl.WdlType.parse_type("Int?"), prefix + "runtime_seconds"),
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
) -> str:
    if not (inp is not None and isinstance(inp, InputSelector)):
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
            base = f'basename({tinp.id()}, "{intype.extension}")'
        else:
            base = f"basename({tinp.id()})"
    else:
        base = tinp.id()

    if intype.optional:
        base = f'if defined({tinp.id()}) then {base} else "generated"'

    return base
