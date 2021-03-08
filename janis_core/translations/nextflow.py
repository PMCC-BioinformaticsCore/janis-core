import os
import json
from typing import Tuple, Dict, List, Optional, Union

from janis_core.types import DataType, Array, String, File, Int, Directory, Stdout, Stderr, Filename, InputSelector, WildcardSelector, Boolean
from janis_core.operators import Operator, StringFormatter, Selector

from janis_core.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument, Tool, ToolType
from janis_core.translations.translationbase import TranslatorBase
from janis_core import Logger
from janis_core.workflow.workflow import StepNode, InputNode, OutputNode, WorkflowBase

import janis_core.translations.nfgen as nfgen


class NextflowTranslator(TranslatorBase):
    LIB_FILENAME = "lib.nf"
    OUTPUT_METADATA_FILENAME = "janis.outputs.metadata"
    NO_FILE_PATH_PREFIX = f"JANIS_NO_FILE"
    PARAM_VAR = "%PARAM%"

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
        # inputsdict = workflow.inputs_map()
        # toolinputs_dict = {k: ToolInput(k, v.intype) for k, v in inputsdict.items()}
        #
        # step_keys = list(workflow.step_nodes.keys())
        # tool_scripts = {}
        # tools = {}
        #
        # tool_scripts[cls.LIB_FILENAME] = cls.generate_generic_functions()
        # for step_id in workflow.step_nodes:
        #     tool = workflow.step_nodes[step_id].tool
        #     tools[tool.versioned_id()] = tool
        #     input_file_prefix = f"{tool.versioned_id()}.input"
        #     # tool_scripts[input_file_prefix] = cls.generate_wf_step_input_vars(tool, step_keys)
        #     tool_scripts[tool.versioned_id()] = cls.tool_script(tool, step_id, input_file_prefix)
        #
        # return cls.workflow_script(workflow, tools), tool_scripts

        imports = [
            cls.init_helper_functions_import()
        ]
        nf_file = nfgen.NFFile(imports=imports, items=[w])

        return nf_file.get_string()

    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> nfgen.process:

        process = cls.init_process(tool)

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

        imports = [
            cls.init_helper_functions_import()
        ]

        items = [
            process,
            cls.init_output_process(tool=tool, nf_process=process),
            cls.init_execution(nf_process=process)
        ]
        nf_file = nfgen.NFFile(imports=imports, items=items)

        return nf_file.get_string()

    @classmethod
    def translate_helper_files(cls, tool):
        helpers = {}

        lib_file = nfgen.NFFile(imports=[], items=cls.generate_generic_functions())

        helpers[cls.LIB_FILENAME] = lib_file.get_string()
        helpers[nfgen.CONFIG_FILENAME] = cls.generate_config()

        return helpers

    @classmethod
    def init_output_process(cls, tool, nf_process: nfgen.NFBase):
        inputs = []
        for o in nf_process.outputs:
            # Always use 'val' qualifier
            inp = nfgen.ProcessInput(qualifier=nfgen.InputProcessQualifier.val, name=nf_process.name + o.name)
            inputs.append(inp)

        tool_outputs = cls.prepare_tool_output(tool)

        script = ""
        for key, val in tool_outputs.items():
            script += f"echo {key}={val} >> {cls.OUTPUT_METADATA_FILENAME}"

        outputs = [nfgen.ProcessOutput(
            qualifier=nfgen.OutputProcessQualifier.path,
            name="janis_output_metadata",
            expression=f"'{cls.OUTPUT_METADATA_FILENAME}'",
        )]

        process = nfgen.Process(
            name="outputs",
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            inputs=inputs,
            outputs=outputs
        )

        return process

    @classmethod
    def init_helper_functions_import(cls):
        items = [nfgen.ImportItem(name=f.name) for f in cls.generate_generic_functions()]
        return nfgen.Import(items, os.path.join(".", cls.DIR_TOOLS, cls.LIB_FILENAME))

    @classmethod
    def init_process(cls, tool):
        # construct script
        script = cls.prepare_script_for_tool(tool)
        pre_script = cls.prepare_expression_inputs(tool)
        pre_script += cls.prepare_input_vars(tool)
        pre_script += cls.prepare_resources_var(tool)

        process = nfgen.Process(
            name=tool.id(),
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            pre_script=pre_script,
        )

        inputs: List[ToolInput] = tool.inputs()
        outputs: List[ToolOutput] = tool.outputs()

        for i in inputs:
            qual = get_input_qualifier_for_inptype(i.input_type)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())

            if isinstance(i.input_type, File):
                inp.as_process_param = f"Channel.fromPath({cls.PARAM_VAR}).collect()"

            process.inputs.append(inp)

        for o in outputs:
            Logger.debug(o.id())
            qual = get_output_qualifier_for_outtype(o.output_type)
            expression = cls.unwrap_expression(o.selector, inputs_dict=tool.inputs_map(), for_output=True)

            # #TODO: make this tidier
            if isinstance(o.output_type, Array):
                if isinstance(o.output_type.subtype(), (File, Directory)):
                    sub_qual = nfgen.OutputProcessQualifier.path
                else:
                    sub_qual = nfgen.OutputProcessQualifier.val

                tuple_elements = expression.strip("][").split(",")
                formatted_list = []
                for expression in tuple_elements:
                    sub_exp = nfgen.TupleElementForOutput(qualifier=sub_qual, expression=expression)
                    formatted_list.append(sub_exp.get_string())

                expression = ", ".join(formatted_list)

            out = nfgen.ProcessOutput(
                qualifier=qual,
                name=o.id(),
                expression=expression,
                is_optional=o.output_type.optional
            )
            process.outputs.append(out)

        return process

    @classmethod
    def init_execution(cls, nf_process: nfgen.NFBase):
        args_list = []
        for i in nf_process.inputs:

            p = f"params.{i.name}"
            if i.as_process_param:
                p = i.as_process_param.replace(cls.PARAM_VAR, p)

            args_list.append(p)

        args = ", ".join(args_list)
        outputs_args = ", ".join(f"{nf_process.name}.out.{o.name}" for o in nf_process.outputs)

        main = [
            f"{nf_process.name}({args})",
            f"outputs({outputs_args})"
        ]

        return nfgen.Workflow(name="", main=main)

    @classmethod
    def workflow_script(cls, workflow: WorkflowBase, tools: Dict[str, Tool]):
        w = nfgen.workflow()

    @classmethod
    def generate_config(cls):
        return f"""
docker.enabled = true
"""

    @classmethod
    def generate_generic_functions(cls):

        functions = [
            nfgen.Function(
                name="optional",
                parameters=["var", "prefix"],
                definition=f"""
var = var.toString()
if (var && ( var != 'None' ) && (! var.contains('{cls.NO_FILE_PATH_PREFIX}')))
{{
    return prefix.toString() + var
}}
else
{{
    return ''
}}
"""
            ),
            nfgen.Function(
                name="boolean_flag",
                parameters=["var", "prefix"],
                definition=f"""
var = var.toString()
if (var == 'True')
{{
    return prefix.toString()
}}
else
{{
    return ''
}}
"""
            ),
        ]

        return functions


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
    def unwrap_expression(
            cls,
            value,
            quote_string=True,
            tool=None,
            for_output=False,
            inputs_dict=None,
            skip_inputs_lookup=False,
            in_shell_script=False,
            **debugkwargs,
    ):
        if value is None:
            if quote_string:
                return "null"
            return None

        if isinstance(value, StepNode):
            raise Exception(
                f"The Step node '{value.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )

        if isinstance(value, list):
            toolid = debugkwargs.get("tool_id", "unwrap_list_expression")
            elements = []
            for i in range(len(value)):
                el = cls.unwrap_expression(
                    value[i],
                    quote_string=quote_string,
                    tool=tool,
                    tool_id=toolid + "." + str(i),
                    inputs_dict=inputs_dict,
                    skip_inputs_lookup=skip_inputs_lookup,
                    for_output=for_output
                )

                elements.append(el)

            list_representation = f"[{', '.join(elements)}]"

            return list_representation

        if isinstance(value, str):
            if quote_string:
                return f"'{value}'"
            else:
                return value
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, Filename):
            formatted = value.generated_filename()
            return formatted

        elif isinstance(value, StringFormatter):
            return cls.translate_string_formatter(
                value,
                in_shell_script=in_shell_script,
                tool=tool,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                **debugkwargs,
            )
        elif isinstance(value, InputSelector):
            if for_output:
                el = cls.prepare_filename_replacements_for(value, inputsdict=inputs_dict)
                return el
            return cls.translate_input_selector(
                selector=value,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                in_shell_script=in_shell_script
            )
        elif isinstance(value, WildcardSelector):
            raise Exception(
                f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
            )
        elif isinstance(value, Operator):
            unwrap_expression_wrap = lambda exp: cls.unwrap_expression(
                exp,
                quote_string=quote_string,
                tool=tool,
                for_output=for_output,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                **debugkwargs,
            )

            return value.to_nextflow(unwrap_expression_wrap, *value.args)

        elif callable(getattr(value, "nextflow", None)):
            return value.nextflow()

        raise Exception(
            "Could not detect type %s to convert to input value" % type(value)
        )

    @classmethod
    def translate_string_formatter(
            cls,
            selector: StringFormatter,
            tool,
            in_shell_script=False,
            inputs_dict=None,
            skip_inputs_lookup=False,
            **debugkwargs,
    ):
        if len(selector.kwargs) == 0:
            return str(selector)

        kwargreplacements = {
            k: f"{cls.unwrap_expression(v, tool=tool, inputs_dict=inputs_dict, skip_inputs_lookup=skip_inputs_lookup, **debugkwargs)}"
            for k, v in selector.kwargs.items()
        }

        arg_val = selector._format
        for k in selector.kwargs:
            arg_val = arg_val.replace(f"{{{k}}}", f"${{{str(kwargreplacements[k])}}}")

        if in_shell_script:
            arg_val = arg_val.replace("\\", "\\\\")

        return arg_val

    @classmethod
    def prepare_filename_replacements_for(cls,
                                          inp: Optional[Selector], inputsdict: Optional[Dict[str, ToolInput]]
                                          ) -> Optional[str]:
        if inp is None or not isinstance(inp, InputSelector):
            return None

        if not inputsdict:
            return f"${inp.input_to_select}.name"
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
                    base = f"{tinp.id()}.simpleName"
                elif tinp.localise_file:
                    base = f"{tinp.id()}.name"
                else:
                    base = f"{tinp.id()}"
            elif isinstance(intype, Filename):
                base = str(cls.unwrap_expression(intype.generated_filename(), inputs_dict=inputsdict, for_output=True))
            elif (
                    intype.is_array()
                    and isinstance(intype.fundamental_type(), (File, Directory))
                    and tinp.localise_file
            ):
                base = f"{tinp.id()}.map(function(el) {{ return el.basename; }})"
            else:
                base = f"{tinp.id()}"

            if intype.optional:
                replacement = f"{tinp.id()} ? {base} : 'generated'"
            else:
                replacement = f"{base}"

            # return f"\"${{{replacement}}}\""
            return replacement

    @classmethod
    def translate_input_selector(cls,
                                 selector: InputSelector,
                                 inputs_dict,
                                 skip_inputs_lookup=False,
                                 in_shell_script=False,
                                 for_output=False
                                 ):
        # TODO: Consider grabbing "path" of File

        sel: str = selector.input_to_select
        if not sel:
            raise Exception("No input was selected for input selector: " + str(selector))

        skip_lookup = skip_inputs_lookup or sel.startswith("runtime_")

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
                        sel = f"{sel}.simpleName"

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

        if in_shell_script:
            sel = f"${sel}"

        return sel

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
        values_provided_from_tool = {
            i.id(): i.default
            for i in tool.tool_inputs()
        }

        count = 0
        inp = {}
        Logger.debug("values_provided_from_tool")
        Logger.debug(values_provided_from_tool)
        for i in tool.tool_inputs():
            val = ad.get(i.id(), values_provided_from_tool.get(i.id()))

            inputsdict = tool.inputs_map()

            if isinstance(i.intype, Filename):
                val = cls.unwrap_expression(i.intype.generated_filename(), inputs_dict=inputsdict)
            elif isinstance(i.intype, Boolean):
                val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""
                if val == "True":
                    val = True
                if val == "False":
                    val = False

            if val is None:
                if isinstance(i.intype, (File, Directory)) \
                        or (isinstance(i.intype, (Array)) and isinstance(i.intype.subtype(), (File, Directory))):
                    count += 1
                    val = f"/{cls.NO_FILE_PATH_PREFIX}{count}"
                else:
                    val = ''

            inp[i.id()] = val

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
        formatted = {}
        for key in inputs:
            # We want list to be formatted as ["xxx", "yyy"] instead of "['xxx', 'yyy']"
            if inputs[key] is not None:
                if type(inputs[key]) is list or type(inputs[key]) is int or type(inputs[key]) is float:
                    val = inputs[key]
                else:
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

        args = [a for a in tool.arguments() or [] if a.position is not None or a.prefix is not None]
        args += [a for a in tool.inputs() or [] if a.position is not None or a.prefix is not None]
        args = sorted(args, key=lambda a: a.position or 0)

        prefix = "  "
        for a in args:
            if isinstance(a, ToolInput):
                pargs.append(f"${a.id()}WithPrefix")
            elif isinstance(a, ToolArgument):
                inputsdict = tool.inputs_map()
                expression = cls.unwrap_expression(
                    a.value,
                    inputsdict=inputsdict,
                    skip_inputs_lookup=True,
                    quote_string=False,
                    in_shell_script=True
                )

                if a.prefix is not None:
                    space = ''
                    if a.separate_value_from_prefix is not False:
                        space = ' '

                    cmd_arg = f"{a.prefix}{space}\"{expression}\""
                else:
                    cmd_arg = expression

                pargs.append(cmd_arg)
            else:
                raise Exception("unknown input type")

        main_script = " \\\n".join(pargs)

        return f"""
echo "{main_script}" > {nfgen.Process.TOOL_EXECUTED_COMMAND_FILENAME}

{main_script} > {nfgen.Process.TOOL_STDOUT_FILENAME}

"""

    @classmethod
    def prepare_tool_output(cls, tool: Tool):
        inputsdict = tool.inputs_map()

        outputs = {}
        for out in tool.outputs():
            if isinstance(out.output_type, Stdout):
                val = "STDOUT"
            elif isinstance(out.output_type, Stderr):
                val = "STDERR"
            else:
                val = f"${tool.id()}{out.tag}"

            outputs[out.tag] = val

        return outputs

    @classmethod
    def prepare_expression_inputs(cls, tool):

        inputsdict = tool.inputs_map()

        script_lines = []
        for i in tool.tool_inputs():
            if isinstance(i.intype, Filename):
                val = cls.unwrap_expression(i.intype.generated_filename(), inputs_dict=inputsdict)

                code = f"""
def {i.id()} = {val}
"""
                script_lines.append(code)

        return "".join(script_lines)

    @classmethod
    def prepare_input_vars(self, tool):
        pre_script_lines = []
        for a in tool.inputs():
            arg_name = ""
            if isinstance(a, ToolInput):
                arg_name = a.id()
            elif isinstance(a, ToolArgument):
                continue
            else:
                raise Exception("unknown input type")

            if isinstance(a.input_type, Array):
                arg_value = f"{arg_name}.join(' ')"
            elif isinstance(a.input_type, (File)):
                arg_value = f"params.{arg_name}"
            else:
                arg_value = arg_name

            prefix = ''
            if a.prefix is not None:
                prefix = a.prefix

            space = ''
            if a.separate_value_from_prefix is not False:
                space = ' '

            if isinstance(a.input_type, Boolean):
                arg = f"boolean_flag({arg_value}, '{prefix}{space}')"
            else:
                if a.input_type.optional:
                    arg = f"optional({arg_value}, '{prefix}{space}')"
                else:
                    arg = f"'{prefix}{space}' + {arg_value}"

            code = f"""
def {arg_name}WithPrefix =  {arg}
"""

            pre_script_lines.append(code)

        return "".join(pre_script_lines)

    @classmethod
    def prepare_resources_var(cls, tool):
        pre_script_lines = []
        for k, v in cls.build_resources_input(
                tool,
                hints=None,
        ).items():
            code = f"""
def {k} =  params.{k}
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
        return nfgen.OutputProcessQualifier.tuple

    elif isinstance(out_type, Stdout):
        return nfgen.OutputProcessQualifier.stdout

    elif isinstance(out_type, (File, Directory)):
        return nfgen.OutputProcessQualifier.path

    return nfgen.OutputProcessQualifier.val
