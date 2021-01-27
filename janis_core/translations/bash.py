import json
from typing import Dict, Tuple, List

from janis_core import Logger
from janis_core.tool.commandtool import ToolArgument, ToolInput, Tool
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

        # Logger.debug(inputsdict)
        # Logger.debug(toolinputs_dict)

        Logger.debug(workflow.get_tools().items())

        tool_commands = {}
        for tool_id, tool in workflow.get_tools().items():
            tool_commands[tool_id] = cls.generate_tool_command(tool)

        Logger.debug(tool_commands)
        command = ""
        for tool_id in tool_commands:
            Logger.debug(tool_id)
            command += f"echo {tool_commands[tool_id]}\n"
            command += f"{command} > $DIR/{tool_id}_stdout 2> $DIR/{tool_id}_stderr\n\n"

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
outputs="${{outputs/STDOUT/$STDOUTPATH}}"
outputs="${{outputs/STDERR/$STDERRPATH}}"
outputs="${{outputs/DIR/$DIR}}"
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

        outputs = {}
        for out in tool.outputs():
            if isinstance(out.output_type, Stdout):
                outputs[out.tag] = "STDOUT"
            elif isinstance(out.output_type, Stderr):
                outputs[out.tag] = "STDERR"
            elif isinstance(out.output_type, File):
                var_name = None
                if isinstance(out.selector, InputSelector):
                    var_name = out.selector.input_to_select
                elif isinstance(out.glob, InputSelector):
                    var_name = out.glob.input_to_select
                elif isinstance(out.selector, Operator):
                    var_name = "..."

                outputs[out.tag] = f"DIR/${var_name}"
            else:
                outputs[out.tag] = "XXX"

        return outputs

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
        #TODO: workflow input
        if tool.type() == ToolType.Workflow:
            for i in tool.tool_inputs():
                return {}


        elif tool.type() == ToolType.CommandTool:
            inp = {}
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

                    val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""

                    if not val:
                        val = []

                    if not isinstance(val, list):
                        val = [val]

                    inp[i.tag] = " ".join(v for v in val)

                    if len(val) > 0:
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
            val = inputs[key]
            lines.append(f"{key}='{val}'")

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
    def unwrap_expression(cls, expression):
        pass

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

    if tool_arg.shell_quote is not False:
        return f"'${name}WithPrefix'"
    else:
        return f"${name}WithPrefix"


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_input.position is not None or tool_input.prefix):
        return None

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
