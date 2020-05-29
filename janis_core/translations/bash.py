from typing import Dict, Tuple, List

from janis_core import Logger
from janis_core.tool.commandtool import ToolArgument, ToolInput
from janis_core.translations import TranslatorBase

from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    StringFormatter,
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
        raise Exception("Not supported for bash translation")

    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
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
        command = " \\\n".join([str_bc, *["  " + a for a in output_args]])

        doc = f"# {tool.id()} bash wrapper"
        meta = tool.bind_metadata() or tool.metadata
        if params_to_include:
            doc += "\n\tNB: this wrapper only contains a subset of the available parameters"
        if meta and meta.documentation:
            doc += "".join(
                "\n# " + l for l in meta.documentation.splitlines(keepends=False)
            )

        return f"""
#!/usr/bin/env sh

{doc}

{command}
        """

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
        workflow,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
    ) -> Dict[str, any]:
        return {}

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf

    @staticmethod
    def stringify_translated_tool(tool):
        return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        return str(inputs)

    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".sh"

    @staticmethod
    def tool_filename(tool):
        return tool.versioned_id() + ".sh"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + ".json"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return None


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
    if tool_arg.shell_quote is not False:
        return f"{tprefix}'${name}'" if tprefix else f"'${name}'"
    else:
        return f"{tprefix}${name}" if tprefix else f"${name}"


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_input.position is not None or tool_input.prefix):
        return None

    name = tool_input.id()
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

    if isinstance(intype, Boolean):
        if tool_input.prefix:
            return tool_input.prefix
        return ""
    elif isinstance(intype, Array):
        Logger.critical("Can't bind arrays onto bash yet")
        return ""

        # expr = name
        #
        # separator = tool_input.separator if tool_input.separator is not None else " "
        # should_quote = isinstance(intype.subtype(), (String, File, Directory))
        # condition_for_binding = None
        #
        # if intype.optional:
        #     expr = f"select_first([{expr}, []])"
        #     condition_for_binding = (
        #         f"(defined({name}) && length(select_first([{name}, []])) > 0)"
        #     )
        #
        # if intype.subtype().optional:
        #     expr = f"select_all({expr})"
        #
        # if should_quote:
        #     if tool_input.prefix_applies_to_all_elements:
        #         separator = f"'{separator}{tprefix} '"
        #     else:
        #         separator = f"'{separator}'"
        #
        #     if tprefix:
        #         expr = f'"{tprefix}\'" + sep("{separator}", {expr}) + "\'"'
        #     else:
        #         expr = f'"\'" + sep("{separator}", {expr}) + "\'"'
        #
        # else:
        #     if tprefix:
        #         expr = f'"{tprefix}" + sep("{separator}", {expr})'
        #     else:
        #         expr = f'sep("{separator}", {expr})'
        # if condition_for_binding is not None:
        #     name = f'~{{if {condition_for_binding} then {expr} else ""}}'
        # else:
        #     name = f"~{{{expr}}}"
    elif (
        isinstance(intype, (String, File, Directory))
        and tool_input.shell_quote is not False
    ):
        return f"{tprefix}'${name}'" if tprefix else f"'${name}'"
        # if tprefix:
        #     # if optional:
        #     # else:
        #     name = f"{tprefix}'${name}'"
        # else:
        #     # if not optional:
        #     # else:
        #     name = f"'${name}'"

    else:
        return f"{tprefix}${name}" if tprefix else f"${name}"
        # if prefix:
        #     if optional:
        #         name = f"~{{if defined({name}) then (\"{tprefix}\" + {name}) else ''}}"
        #     else:
        #         name = f"{tprefix}~{{{name}}}"
        # else:
        #     name = f"~{{{name}}}"


if __name__ == "__main__":
    from janis_unix.tools import Echo

    Echo().translate("shell")
