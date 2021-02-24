import json
from typing import Tuple, Dict, List, Optional, Union

from janis_core.types import DataType, Array, String, File, Int, Directory, Stdout, Stderr

from janis_core.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument, Tool, ToolType
from janis_core.translations.translationbase import TranslatorBase
from janis_core import Logger

import janis_core.translations.nfgen as nfgen


class NextflowTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="nextflow")

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
        pre_script = cls.prepare_optional_inputs(tool)

        process = nfgen.Process(
            name=tool.id(),
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            pre_script=pre_script
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

        return process.get_string()

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        raise Exception("CodeTool is not currently supported in Nextflow translation")

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
        # if tool.type() == ToolType.Workflow:
        #     values_provided_from_tool = {
        #         i.id(): i.value or i.default
        #         for i in tool.input_nodes.values()
        #         if i.value or (i.default and not isinstance(i.default, Selector))
        #     }

        # inp = {
        #     i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
        #     for i in tool.tool_inputs()
        #     if i.default is not None
        #        or not i.intype.optional
        #        or i.id() in ad
        #        or i.id() in values_provided_from_tool
        # }

        # inp = {
        #     i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
        #     for i in
        #     tool.tool_inputs()
        # }

        count = 0
        inp = {}
        for i in tool.tool_inputs():
            val = ad.get(i.id(), values_provided_from_tool.get(i.id()))

            if val is None:
                if isinstance(i.intype, (File, Directory)):
                    count += 1
                    val = f"$baseDir/no_file{count}"

            inp[i.id()] = val


        # if merge_resources:
        #     for k, v in cls.build_resources_input(
        #             tool,
        #             hints,
        #             max_cores=max_cores,
        #             max_mem=max_mem,
        #             max_duration=max_duration,
        #     ).items():
        #         inp[k] = ad.get(k, v)

        return inp

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf

    @staticmethod
    def stringify_translated_tool(tool):
        return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        formatted = {}
        for key in inputs:
            if inputs[key] is not None:
                val = str(inputs[key])
            else:
                val = ''

            formatted[key] = val

        return json.dumps(formatted)

    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".nf"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.versioned_id() + ".input.json"

    @staticmethod
    def tool_filename(tool):
        prefix = tool
        if isinstance(tool, Tool):
            prefix = tool.versioned_id()

        return prefix + ".nf"


    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass

    @classmethod
    def prepare_script_for_tool(cls, tool: CommandTool):
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
            if isinstance(a, ToolInput):
                pargs.append(f"${a.id()}WithPrefix")
            elif isinstance(a, ToolArgument):
                pargs.append(f"{a.value}")
            else:
                raise Exception("unknown input type")

            # if isinstance(a, ToolInput):
            #     pargs.append(prefix + f"{a.prefix or ''} !{{{a.id()}}}")
            # elif isinstance(a, ToolArgument):
            #     pargs.append(prefix + f"{a.prefix or ''} !{{{a.value}}}")

        main_script = " \\\n".join(pargs)
        outputs = cls.prepare_tool_output(tool)
        esc = '\\\\"'

        return f"""
{main_script}

export DIR=\$(pwd)
export STDOUTPATH=\$(pwd)/.command.out
export STDERRPATH=\$(pwd)/.command.err

outputs="{json.dumps(outputs).replace('"', esc)}"
outputs="\${{outputs//STDOUT/\$STDOUTPATH}}"
outputs="\${{outputs//STDERR/\$STDERRPATH}}"
outputs="\${{outputs//DIR/\$DIR}}"
echo \$outputs > .janis.out
"""

    @classmethod
    def prepare_tool_output(cls, tool: Tool):
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
                if sel is None:
                    sel = out.glob

                val = "DIR/" + cls.unwrap_expression(sel, inputs_dict=inputsdict)

            else:
                val = "XXX"

            outputs[out.tag] = val

        return outputs


    @classmethod
    def prepare_optional_inputs(self, tool):
        pre_script_lines = []
        for a in tool.inputs():
            arg_name = ""
            if isinstance(a, ToolInput):
                arg_name = a.id()
            elif isinstance(a, ToolArgument):
                continue
            else:
                raise Exception("unknown input type")

            code = f"""
def {arg_name}WithPrefix =  params.{arg_name} && ( params.{arg_name} != 'None' ) && (! params.{arg_name}.contains('/no_file')) ? "{a.prefix or ''} $params.{arg_name}" : ''
"""

            pre_script_lines.append(code)

        return "".join(pre_script_lines)


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
