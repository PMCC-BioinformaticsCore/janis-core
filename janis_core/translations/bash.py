import json
from typing import Dict, Tuple, List, Optional

from janis_core import Logger
from janis_core.tool.commandtool import ToolArgument, ToolInput, Tool, ToolOutput
from janis_core.workflow.workflow import StepNode, InputNode, OutputNode

from janis_core.tool.tool import ToolType
from janis_core.translations import TranslatorBase

from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    String,
    Selector,
    Directory,
    Stdout,
    Stderr,
    Array,
    Boolean,
    Filename,
    File,
)

from janis_core.operators import Operator

class BashTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="bash")

    @classmethod
    def translate_workflow(
        cls,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:

        # raise Exception("Not supported for bash translation")

        inputsdict = workflow.inputs_map()
        toolinputs_dict = {k: ToolInput(k, v.intype) for k, v in inputsdict.items()}

        tool_commands = {}
        for tool_id, tool in workflow.get_tools().items():
            tool_commands[tool_id] = cls.generate_tool_command(tool)

        command = ""
        for tool_id in tool_commands:
            command += f"echo \"{tool_commands[tool_id]}\"\n"
            command += f"{tool_commands[tool_id]} > $DIR/{tool_id}_stdout 2> $DIR/{tool_id}_stderr\n\n"

        return (f"""
#!/usr/bin/env sh

source $1

DIR=$(pwd)
STDOUTPATH=$(pwd)/stdout
STDERRPATH=$(pwd)/stderr

{command}

# END
""", [])


    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        params_to_include = None
        if tool.connections:
            params_to_include = set(tool.connections.keys())

        doc = f"# {tool.id()} bash wrapper"
        meta = tool.bind_metadata() or tool.metadata
        if params_to_include:
            doc += "\n\tNB: this wrapper only contains a subset of the available parameters"
        if meta and meta.documentation:
            doc += "".join(
                "\n# " + l for l in meta.documentation.splitlines(keepends=False)
            )

        outputs = cls.generate_outputs(tool)
        esc = '\\"'

        command = cls.generate_tool_command(tool)

        return f"""
#!/usr/bin/env sh

source $1

DIR=$(pwd)
STDOUTPATH=$(pwd)/stdout
STDERRPATH=$(pwd)/stderr

{doc}

echo {command}
{command} > $STDOUTPATH 2> $STDERRPATH

outputs="{json.dumps(outputs).replace('"', esc)}"
outputs="${{outputs//STDOUT/$STDOUTPATH}}"
outputs="${{outputs//STDERR/$STDERRPATH}}"
outputs="${{outputs//DIR/$DIR}}"
echo $outputs
"""

    @classmethod
    def generate_tool_command(cls, tool):
        args: List[ToolArgument] = sorted(
            [*(tool.arguments() or []), *(tool.inputs() or [])],
            key=lambda a: (a.position or 0),
        )

        params_to_include = None
        if tool.connections:
            params_to_include = set(tool.connections.keys())

        bc = tool.base_command()
        if bc is None:
            bc = []
        elif not isinstance(bc, list):
            bc = [bc]

        output_args = []
        for a in args:
            if isinstance(a, ToolInput):
                if params_to_include and a.id() not in params_to_include:
                    # skip if we're limiting to specific commands
                    continue
                arg = translate_command_input(tool_input=a, inputsdict={})
                if not arg:
                    Logger.warn(f"Parameter {a.id()} was skipped")
                    continue
                output_args.append(arg)
            else:
                output_args.append(
                    translate_command_argument(tool_arg=a, inputsdict={})
                )

        str_bc = " ".join(f"'{c}'" for c in bc)
        command = " \\\n".join([str_bc, *[a for a in output_args]])

        return command

    @classmethod
    def generate_outputs(cls, tool: Tool):

        inputsdict = tool.inputs_map()

        outputs = {}
        for out in tool.outputs():
            if isinstance(out.output_type, Array):
                val = []
                for sel in out.selector:
                    val.append("DIR/" + cls.unwrap_expression(sel, inputs_dict=inputsdict))

            elif isinstance(out.output_type, Stdout):
                val = "STDOUT"
            elif isinstance(out.output_type, Stderr):
                val = "STDERR"
            elif isinstance(out.output_type, File):
                sel = out.selector
                if not sel:
                    sel = out.glob

                val = "DIR/" + cls.unwrap_expression(sel, inputs_dict=inputsdict)

            else:
                val = "XXX"

            outputs[out.tag] = val

        return outputs

    @classmethod
    def read_selector(cls, sel: Selector):
        var_name = None
        if isinstance(sel, InputSelector):
            var_name = sel.input_to_select
        # elif isinstance(out.glob, InputSelector):
        #     var_name = out.glob.input_to_select
        elif isinstance(sel, Operator):
            var_name = "..."

        return f"DIR/${var_name}"

    @classmethod
    def unwrap_expression(
            cls,
            value,
            code_environment=True,
            selector_override=None,
            tool=None,
            for_output=False,
            inputs_dict=None,
            **debugkwargs,
    ):
        if value is None:
            if code_environment:
                return "null"
            return None

        if isinstance(value, StepNode):
            raise Exception(
                f"The Step node '{value.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )

        if isinstance(value, list):
            toolid = debugkwargs.get("tool_id", "unwrap_list_expression")
            inner = ", ".join(
                cls.unwrap_expression(
                    value[i],
                    code_environment=True,
                    selector_override=selector_override,
                    tool=tool,
                    tool_id=toolid + "." + str(i),
                    inputs_dict=inputs_dict,
                )
                for i in range(len(value))
            )
            return cls.wrap_in_codeblock_if_required(
                f"[{inner}]", is_code_environment=code_environment
            )

        if isinstance(value, str):
            return value
            # if not code_environment:
            #     return value
            # return cls.quote_values_if_code_environment(
            #     cls.prepare_escaped_string(value), code_environment
            # )
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, Filename):
            # value.generated_filenamecwl() if code_environment else f"$({value.generated_filenamecwl()})"
            return cls.quote_values_if_code_environment(
                value.generated_filename(), code_environment
            )
        # elif isinstance(value, AliasSelector):
        #     return cls.unwrap_expression(
        #         value.inner_selector,
        #         code_environment=code_environment,
        #         selector_override=selector_override,
        #         inputs_dict=inputs_dict,
        #         for_output=for_output,
        #         tool=tool,
        #         **debugkwargs,
        #     )
        #
        # elif isinstance(value, StringFormatter):
        #     return translate_string_formatter(
        #         value,
        #         selector_override=selector_override,
        #         code_environment=code_environment,
        #         tool=tool,
        #         inputs_dict=inputs_dict,
        #         **debugkwargs,
        #     )
        # elif isinstance(value, InputNodeSelector):
        #     return translate_input_selector(
        #         InputSelector(value.id()),
        #         code_environment=code_environment,
        #         selector_override=selector_override,
        #         inputs_dict=inputs_dict,
        #         skip_inputs_lookup=True,
        #     )
        # elif isinstance(value, StepOutputSelector):
        #     sel = f"{value.node.id()}/{value.tag}"
        #     if sel in selector_override:
        #         return selector_override[sel]
        #     raise Exception(
        #         "An internal error occurred when unwrapping an operator, found StepOutputSelector with no alias"
        #     )
        # elif isinstance(value, ResourceSelector):
        #     if not tool:
        #         raise Exception(
        #             f"Tool must be provided when unwrapping ResourceSelector: {type(value).__name__}"
        #         )
        #     operation = value.get_operation(tool, hints={})
        #     return cls.unwrap_expression(
        #         operation,
        #         code_environment=code_environment,
        #         tool=tool,
        #         inputs_dict=inputs_dict,
        #         **debugkwargs,
        #     )
        #
        # elif for_output and isinstance(value, (Stderr, Stdout)):
        #     # next few ones we rely on the globs being
        #     if isinstance(value, Stdout):
        #         return "self[0]"
        #     elif isinstance(value, Stderr):
        #         return "self[1]"

        elif isinstance(value, InputSelector):
            if for_output:
                el = cls.prepare_filename_replacements_for(value, inputsdict=inputs_dict)
                return cls.wrap_in_codeblock_if_required(
                    el, is_code_environment=code_environment
                )
            return cls.translate_input_selector(
                selector=value,
                code_environment=code_environment,
                selector_override=selector_override,
                inputs_dict=inputs_dict,
            )
        elif isinstance(value, WildcardSelector):
            raise Exception(
                f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
            )
        elif isinstance(value, Operator):
            unwrap_expression_wrap = lambda exp: cls.unwrap_expression(
                exp,
                code_environment=True,
                selector_override=selector_override,
                tool=tool,
                for_output=for_output,
                inputs_dict=inputs_dict,
                **debugkwargs,
            )
            return cls.wrap_in_codeblock_if_required(
                value.to_shell(unwrap_expression_wrap, *value.args),
                is_code_environment=code_environment,
            )
        elif callable(getattr(value, "cwl", None)):
            return value.cwl()
        # elif isinstance(value, Operator):

        raise Exception(
            "Could not detect type %s to convert to input value" % type(value)
        )

    @classmethod
    def translate_input_selector(cls,
            selector: InputSelector,
            code_environment,
            inputs_dict,
            selector_override=None,
            skip_inputs_lookup=False,
    ):
        # TODO: Consider grabbing "path" of File

        sel: str = selector.input_to_select
        if not sel:
            raise Exception("No input was selected for input selector: " + str(selector))

        skip_lookup = skip_inputs_lookup or sel.startswith("runtime_")

        if selector_override and sel in selector_override:
            sel = selector_override[sel]
        else:
            sel = f"${sel}"

        if not skip_lookup:

            if inputs_dict is None:
                raise Exception(
                    f"An internal error occurred when translating input selector '{sel}': the inputs dictionary was None"
                )
            if selector.input_to_select not in inputs_dict:
                raise Exception(
                    f"Couldn't find the input '{sel}' for the InputSelector(\"{sel}\")"
                )

            tinp: ToolInput = inputs_dict[selector.input_to_select]

            intype = tinp.intype
            if selector.remove_file_extension:
                if intype.is_base_type((File, Directory)):
                    potential_extensions = (
                        intype.get_extensions() if intype.is_base_type(File) else None
                    )
                    if selector.remove_file_extension and potential_extensions:
                        sel = f"{sel}.basename"
                        for ext in potential_extensions:
                            sel += f'.replace(/{ext}$/, "")'

                elif intype.is_array() and isinstance(
                        intype.fundamental_type(), (File, Directory)
                ):
                    inner_type = intype.fundamental_type()
                    extensions = (
                        inner_type.get_extensions()
                        if isinstance(inner_type, File)
                        else None
                    )

                    inner_sel = f"el.basename"
                    if extensions:
                        for ext in extensions:
                            inner_sel += f'.replace(/{ext}$/, "")'
                    sel = f"{sel}.map(function(el) {{ return {inner_sel}; }})"
                else:
                    Logger.warn(
                        f"InputSelector {sel} is requesting to remove_file_extension but it has type {tinp.input_type.id()}"
                    )
            # elif tinp.localise_file:
            #     if intype.is_base_type((File, Directory)):
            #         sel += ".basename"
            #     elif intype.is_array() and isinstance(
            #             intype.fundamental_type(), (File, Directory)
            #     ):
            #         sel = f"{sel}.map(function(el) {{ return el.basename; }})"

        return sel if code_environment else f"$({sel})"

    @classmethod
    def prepare_filename_replacements_for(cls,
            inp: Optional[Selector], inputsdict: Optional[Dict[str, ToolInput]]
    ) -> Optional[str]:
        if inp is None or not isinstance(inp, InputSelector):
            return None

        if not inputsdict:
            return "inputs." + inp.input_to_select + ".basename"
            # raise Exception(
            #     f"Couldn't generate filename as an internal error occurred (inputsdict did not contain {inp.input_to_select})"
            # )

        if isinstance(inp, InputSelector):
            if inp.input_to_select not in inputsdict:
                raise Exception(
                    f"The InputSelector '{inp.input_to_select}' did not select a valid input"
                )

            tinp = inputsdict.get(inp.input_to_select)
            intype = tinp.intype

            if intype.is_base_type((File, Directory)):
                potential_extensions = (
                    intype.get_extensions() if intype.is_base_type(File) else None
                )
                if inp.remove_file_extension and potential_extensions:
                    base = f"inputs.{tinp.id()}.basename"
                    for ext in potential_extensions:
                        base += f'.replace(/{ext}$/, "")'
                elif tinp.localise_file:
                    base = f"inputs.{tinp.id()}.basename"
                else:
                    base = f"inputs.{tinp.id()}"
            elif (
                    intype.is_array()
                    and isinstance(intype.fundamental_type(), (File, Directory))
                    and tinp.localise_file
            ):
                base = f"inputs.{tinp.id()}.map(function(el) {{ return el.basename; }})"
            else:
                base = "inputs." + tinp.id()

            if intype.optional:
                replacement = f'inputs.{tinp.id()} ? {base} : "generated"'
            else:
                replacement = f"{base}"

            return replacement

    @classmethod
    def wrap_in_codeblock_if_required(cls, value, is_code_environment):
        return value if is_code_environment else f"$({value})"

    @classmethod
    def quote_values_if_code_environment(cls, value, is_code_environment):
        return f'"{value}"' if is_code_environment else value

    @classmethod
    def prepare_escaped_string(cls, value: str):
        return json.dumps(value)[1:-1]

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        raise Exception("CodeTool is not currently supported in bash translation")

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

        ad = additional_inputs or {}
        values_provided_from_tool = {}
        if tool.type() == ToolType.Workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value or (i.default and not isinstance(i.default, Selector))
            }
            # TODO: add mapping from inputdict

        elif tool.type() == ToolType.CommandTool:
            values_provided_from_tool = {
                i.id(): i.default
                for i in tool.tool_inputs()
            }

        # inp = {
        #     i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
        #     for i in tool.tool_inputs()
        #     if i.default is not None
        #        or not i.intype.optional
        #        or i.id() in ad
        #        or i.id() in values_provided_from_tool
        # }

        # Build input variables and another copy of each input variable with its prefix attached
        inp = {}
        #TODO: workflow input
        if tool.type() == ToolType.Workflow:
            for i in tool.tool_inputs():
                if i.default is not None \
                        or not i.intype.optional \
                        or i.id() in ad \
                        or i.id() in values_provided_from_tool:
                    val = ad.get(i.id(), values_provided_from_tool.get(i.id()))

                    if not val:
                        val = []

                    if not isinstance(val, list):
                        val = [val]

                    inp[i.tag] = " ".join(str(v) for v in val)
                    inp[i.tag + "WithPrefix"] = inp[i.tag]

        elif tool.type() == ToolType.CommandTool:
            for i in tool.inputs():
                if i.default is not None \
                   or not i.input_type.optional \
                   or i.tag in ad \
                   or i.tag in values_provided_from_tool:

                    prefix = i.prefix if i.prefix else ""
                    tprefix = prefix

                    if prefix and i.separate_value_from_prefix:
                        tprefix += " "

                    ad.get(i.tag)
                    values_provided_from_tool.get(i.tag)

                    inputsdict = tool.inputs_map()

                    if isinstance(i.input_type, Filename):
                        val = cls.unwrap_expression(i.input_type.generated_filename(), inputs_dict=inputsdict)
                    else:
                        val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""

                    if val == "":
                        val = []

                    if not isinstance(val, list):
                        val = [val]

                    inp[i.tag] = " ".join(str(v) for v in val)

                    if len(val) > 0 and (i.prefix or i.position):
                        inp[i.tag + "WithPrefix"] = " ".join(tprefix + v for v in val)
                    else:
                        inp[i.tag + "WithPrefix"] = ""

        if merge_resources:
            for k, v in cls.build_resources_input(
                    tool,
                    hints,
                    max_cores=max_cores,
                    max_mem=max_mem,
                    max_duration=max_duration,
            ).items():
                inp[k] = ad.get(k, v)

        return inp

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf

    @staticmethod
    def stringify_translated_tool(tool):
        return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        lines = []
        Logger.debug(f"inputs: {inputs}")
        for key in inputs:
            val = str(inputs[key])
            lines.append(f"{key}=\"{val}\"")

        return "\n".join(lines)


    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".sh"

    @staticmethod
    def tool_filename(tool):
        return tool.versioned_id() + ".sh"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + ".input.sh"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return None


    @classmethod
    def stdout_output_name(cls, tool):
        stdout_outputs = []
        for o in tool.outputs():
            if isinstance(o.output_type, Stdout):
                stdout_outputs.append(o.tag)

        if not stdout_outputs:
            return None

        if len(stdout_outputs) != 1:
            raise Exception("There is more than out one output with type Stdout")

        return stdout_outputs[0]


def translate_command_argument(tool_arg: ToolArgument, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_arg.position is not None or tool_arg.prefix):
        return None

    separate_value_from_prefix = tool_arg.separate_value_from_prefix is not False
    prefix = tool_arg.prefix if tool_arg.prefix else ""
    tprefix = prefix

    if prefix and separate_value_from_prefix:
        tprefix += " "

    name = tool_arg.value
    # if tool_arg.shell_quote is not False:
    #     return f"{tprefix}'${name}'" if tprefix else f"'${name}'"
    # else:
    #     return f"{tprefix}${name}" if tprefix else f"${name}"
    # if tool_arg.shell_quote is not False:
    #     return f"'${name}'" if tprefix else f"'${name}'"
    # else:
    #     return f"${name}" if tprefix else f"${name}"

    # if tool_arg.shell_quote is not False:
    #     return f"'${name}WithPrefix'"
    # else:
    #     return f"${name}WithPrefix"

    arg_val = f"{name}"
    arg_val = arg_val.replace("inputs.", "$")

    if tool_arg.shell_quote is not False:
        return f"'{arg_val}'"
    else:
        return f"{arg_val}"


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    # if not (tool_input.position is not None or tool_input.prefix):
    #     return None

    name = tool_input.id()
    return f"${name}WithPrefix"

    # name = tool_input.id()
    # intype = tool_input.input_type
    #
    # optional = (not isinstance(intype, Filename) and intype.optional) or (
    #     isinstance(tool_input.default, CpuSelector) and tool_input.default is None
    # )
    # position = tool_input.position
    #
    # separate_value_from_prefix = tool_input.separate_value_from_prefix is not False
    # prefix = tool_input.prefix if tool_input.prefix else ""
    # tprefix = prefix
    #
    # intype = tool_input.input_type
    #
    # is_flag = isinstance(intype, Boolean)
    #
    # if prefix and separate_value_from_prefix and not is_flag:
    #     tprefix += " "
    #
    # if isinstance(intype, Boolean):
    #     # if tool_input.prefix:
    #     #     return tool_input.prefix
    #     # return ""
    #     return f"${name}WithPrefix"
    # elif isinstance(intype, Array):
    #     # Logger.critical("Can't bind arrays onto bash yet")
    #     # return ""
    #
    #     return f"${name}WithPrefix"
    #
    #     # expr = name
    #     #
    #     # separator = tool_input.separator if tool_input.separator is not None else " "
    #     # should_quote = isinstance(intype.subtype(), (String, File, Directory))
    #     # condition_for_binding = None
    #     #
    #     # if intype.optional:
    #     #     expr = f"select_first([{expr}, []])"
    #     #     condition_for_binding = (
    #     #         f"(defined({name}) && length(select_first([{name}, []])) > 0)"
    #     #     )
    #     #
    #     # if intype.subtype().optional:
    #     #     expr = f"select_all({expr})"
    #     #
    #     # if should_quote:
    #     #     if tool_input.prefix_applies_to_all_elements:
    #     #         separator = f"'{separator}{tprefix} '"
    #     #     else:
    #     #         separator = f"'{separator}'"
    #     #
    #     #     if tprefix:
    #     #         expr = f'"{tprefix}\'" + sep("{separator}", {expr}) + "\'"'
    #     #     else:
    #     #         expr = f'"\'" + sep("{separator}", {expr}) + "\'"'
    #     #
    #     # else:
    #     #     if tprefix:
    #     #         expr = f'"{tprefix}" + sep("{separator}", {expr})'
    #     #     else:
    #     #         expr = f'sep("{separator}", {expr})'
    #     # if condition_for_binding is not None:
    #     #     name = f'~{{if {condition_for_binding} then {expr} else ""}}'
    #     # else:
    #     #     name = f"~{{{expr}}}"
    # elif (
    #     isinstance(intype, (String, File, Directory))
    #     and tool_input.shell_quote is not False
    # ):
    #     # return f"{tprefix}\"${name}\"" if tprefix else f"\"${name}\""
    #     return f"'${name}WithPrefix'"
    #     # if tprefix:
    #     #     # if optional:
    #     #     # else:
    #     #     name = f"{tprefix}'${name}'"
    #     # else:
    #     #     # if not optional:
    #     #     # else:
    #     #     name = f"'${name}'"
    #
    # else:
    #     # return f"{tprefix}${name}" if tprefix else f"${name}"
    #     return f"${name}WithPrefix"
    #     # if prefix:
    #     #     if optional:
    #     #         name = f"~{{if defined({name}) then (\"{tprefix}\" + {name}) else ''}}"
    #     #     else:
    #     #         name = f"{tprefix}~{{{name}}}"
    #     # else:
    #     #     name = f"~{{{name}}}"


if __name__ == "__main__":
    from janis_unix.tools import Echo

    Echo().translate("shell")
