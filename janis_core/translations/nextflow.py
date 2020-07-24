from typing import Tuple, Dict, List, Optional, Union

from janis_core.types import DataType, Array, String, File, Int, Directory, Stdout

from janis_core.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument
from janis_core.translations.translationbase import TranslatorBase

import janis_core.translations.nfgen as nfgen


class NextflowTranslator(TranslatorBase):
    @classmethod
    def translate_workflow(
        cls,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:
        pass

    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> nfgen.process:
        # construct script
        script = cls.prepare_script_for_tool(tool)
        process = nfgen.Process(
            name=tool.id(),
            script=script,
            script_type=nfgen.ProcessScriptType.shell,
            script_quote="'",
        )

        inputs: List[ToolInput] = tool.inputs()
        outputs: List[ToolOutput] = tool.outputs()
        inpmap = {i.id(): i for i in inputs}

        for i in inputs:
            qual = get_input_qualifier_for_inptype(i.input_type)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())
            process.inputs.append(inp)

        for o in outputs:
            qual = get_output_qualifier_for_outtype(o.output_type)
            out = nfgen.ProcessOutput(
                qualifier=qual, name=o.id(), is_optional=o.output_type.optional
            )
            process.outputs.append(out)

        if with_container:
            container = (
                NextflowTranslator.get_container_override_for_tool(
                    tool, container_override
                )
                or tool.container()
            )

            if container is not None:
                process.directives.append(
                    nfgen.ContainerDirective(
                        cls.unwrap_expression(container, is_code_environment=True)
                    )
                )
            elif not allow_empty_container:
                raise Exception(
                    f"The tool '{tool.id()}' did not have a container and no container override was specified. "
                    f"Although not recommended, Janis can export empty docker containers with the parameter "
                    f"'allow_empty_container=True' or --allow-empty-container"
                )

        return process

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        pass

    @classmethod
    def prepare_string_if_required(cls, value, is_code_environment):
        return f'"{value}"' if is_code_environment else value

    @classmethod
    def unwrap_expression(cls, expression, is_code_environment=True):
        if isinstance(expression, str):
            return cls.prepare_string_if_required(expression, is_code_environment)

        raise Exception(
            f"Could not detect type '{type(expression)}' to unwrap to nextflow"
        )

    @staticmethod
    def stringify_translated_workflow(wf):
        pass

    @staticmethod
    def stringify_translated_tool(tool):
        pass

    @staticmethod
    def stringify_translated_inputs(inputs):
        pass

    @staticmethod
    def workflow_filename(workflow):
        pass

    @staticmethod
    def inputs_filename(workflow):
        pass

    @staticmethod
    def tool_filename(tool):
        pass

    @staticmethod
    def resources_filename(workflow):
        pass

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass

    @staticmethod
    def prepare_script_for_tool(tool: CommandTool):
        bc = tool.base_command()
        pargs = []

        if bc:
            pargs.append(" ".join(bc) if isinstance(bc, list) else str(bc))

        args = sorted(
            [*(tool.arguments() or []), *(tool.inputs() or [])],
            key=lambda a: a.position,
        )

        prefix = "  "
        for a in args:
            # if isinstance(a, ToolInput):
            pargs.append(prefix + f"{a.prefix or ''} !{{{a.id()}}}")

        return " \\\n".join(pargs)


def get_input_qualifier_for_inptype(inp_type: DataType) -> nfgen.InputProcessQualifier:

    if isinstance(inp_type, Array):
        inp_type = inp_type.fundamental_type()

    if isinstance(inp_type, (File, Directory)):
        return nfgen.InputProcessQualifier.path
    return nfgen.InputProcessQualifier.val


def get_output_qualifier_for_outtype(
    out_type: DataType,
) -> nfgen.OutputProcessQualifier:
    if isinstance(out_type, Array):
        out_type = out_type.fundamental_type()

    if isinstance(out_type, Stdout):
        return nfgen.OutputProcessQualifier.stdout

    if isinstance(out_type, (File, Directory)):
        return nfgen.OutputProcessQualifier.path
    return nfgen.OutputProcessQualifier.val
