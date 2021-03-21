import os
import posixpath
import json
from typing import Tuple, Dict, List, Optional, Union, Any

from janis_core.types import DataType, Array, String, File, Int, Float, Double, Directory, Stdout, Stderr, Filename, InputSelector, WildcardSelector, Boolean
from janis_core.operators import Operator, StringFormatter, Selector

from janis_core.tool.commandtool import CommandTool, ToolInput, ToolOutput, ToolArgument, Tool, ToolType, TOutput, TInput
from janis_core.code.codetool import CodeTool
from janis_core.code.pythontool import PythonTool
from janis_core.translations.translationbase import TranslatorBase
from janis_core import Logger
from janis_core.workflow.workflow import StepNode, InputNode, OutputNode, WorkflowBase
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core.translationdeps.supportedtranslations import SupportedTranslation


import janis_core.translations.nfgen as nfgen


class NextflowTranslator(TranslatorBase):
    LIB_FILENAME = "lib.nf"
    OUTPUT_METADATA_FILENAME = "janis.outputs.metadata"
    NO_FILE_PATH_PREFIX = f"JANIS_NO_FILE"
    PARAM_VAR = "%PARAM%"
    CHANNEL_PARAM = "%CHANNEL_PARAM%"
    PYTHON_CODE_FILE_PATH_PARAM = "%PYTHON_CODE_FILE_PATH%"
    PYTHON_CODE_OUTPUT_FILENAME_PREFIX = "janis_out_"

    INPUT_IN_SELECTORS = {}

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
        inputsdict = workflow.inputs_map()
        toolinputs_dict = {k: ToolInput(k, v.intype) for k, v in inputsdict.items()}

        step_keys = list(workflow.step_nodes.keys())
        nf_items = {}
        process_files = {}
        tools = {}

        wf_outputs = cls.generate_wf_tool_outputs(workflow)
        # Logger.debug("wf_outputs")
        # Logger.debug(wf_outputs)

        for step_id in workflow.step_nodes:
            tool = workflow.step_nodes[step_id].tool
            tools[tool.versioned_id()] = tool

            provided_inputs = cls.generate_wf_tool_inputs(tool, step_keys)
            # Logger.debug(f"{step_id}: provided_inputs")
            # Logger.debug(provided_inputs)

            if isinstance(tool, CommandTool):
                nf_item = cls.init_process_command_tool(tool, step_id, provided_inputs)
                nf_item = cls.handle_container(tool, nf_item, with_container, allow_empty_container, container_override)
            if isinstance(tool, PythonTool):
                nf_item = cls.init_process_python_code_tool(tool, step_id, provided_inputs)
                nf_item = cls.handle_container(tool, nf_item, with_container, allow_empty_container, container_override)
            elif isinstance(tool, WorkflowBase):
                nf_item = cls.init_subworkflow(tool, step_id, provided_inputs)
            elif isinstance(tool, CodeTool):
                raise Exception("Only PythonTool code tool is supported for the moment")

            nf_items[tool.versioned_id()] = nf_item

            nf_file = nfgen.NFFile(imports=[cls.init_helper_functions_import(".")], items=[nf_item])
            process_files[tool.versioned_id()] = nf_file

        out_process_inp, out_process_out = cls.prepare_output_process_params_for_worfklow(workflow)

        imports = [cls.init_helper_functions_import()] + cls.init_tool_steps_import(nf_items)
        items = [
            cls.init_output_process(out_process_inp, out_process_out),
            cls.init_workflow(workflow=workflow, nf_items=nf_items)
        ]
        nf_file = nfgen.NFFile(imports=imports, items=items)

        tool_scripts = {t: process_files[t].get_string() for t in process_files}

        return (nf_file.get_string(), tool_scripts)

    @classmethod
    def init_subworkflow(cls, worfklow: WorkflowBase, step_id: str, provided_inputs):
        pass

    @classmethod
    def handle_container(
            cls,
            tool,
            process: nfgen.Process,
            with_container: bool = True,
            allow_empty_container=False,
            container_override: dict = None
    ):
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
                    cls.unwrap_expression(container, is_code_environment=True, tool=tool)
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
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> nfgen.process:

        process = cls.init_process_command_tool(tool)
        process = cls.handle_container(tool, process, with_container, allow_empty_container, container_override)

        imports = [
            cls.init_helper_functions_import()
        ]

        out_process_inp, out_process_out = cls.prepare_output_process_params_for_tool(tool, process)

        items = [
            process,
            cls.init_output_process(inputs=out_process_inp, tool_outputs=out_process_out),
            cls.init_tool_execution(tool, nf_process=process)
        ]
        nf_file = nfgen.NFFile(imports=imports, items=items)

        return nf_file.get_string()

    @classmethod
    def translate_code_tool_internal(
            cls,
            tool,
            with_container=True,
            allow_empty_container=False,
            container_override: dict = None,
    ):
        # raise Exception("CodeTool is not currently supported in Nextflow translation")

        if isinstance(tool, PythonTool):
            process = cls.init_process_python_code_tool(tool)
            process = cls.handle_container(tool, process, with_container, allow_empty_container, container_override)

            imports = [
                cls.init_helper_functions_import()
            ]

            out_process_inp, out_process_out = cls.prepare_output_process_params_for_tool(tool, process)

            items = [
                process,
                cls.init_output_process(inputs=out_process_inp, tool_outputs=out_process_out),
                cls.init_tool_execution(tool, nf_process=process)
            ]
            nf_file = nfgen.NFFile(imports=imports, items=items)

            return nf_file.get_string()
        else:
            raise Exception("Only PythonTool code tool is supported for the moment.")

    @classmethod
    def generate_python_script(cls, tool: PythonTool):
        return tool.prepared_script(SupportedTranslation.Nextflow)

    @classmethod
    def translate_helper_files(cls, tool):
        helpers = {}

        lib_file = nfgen.NFFile(imports=[], items=cls.generate_generic_functions())

        helpers[cls.LIB_FILENAME] = lib_file.get_string()
        helpers[nfgen.CONFIG_FILENAME] = cls.generate_config()

        # NOTE: for python tool inside a workflow, it will be handled in translate_workflow
        if isinstance(tool, PythonTool):
            helpers["__init__.py"] = ""
            helpers[f"{tool.versioned_id()}.py"] = cls.generate_python_script(tool)

        return helpers

    @classmethod
    def prepare_output_process_params_for_tool(cls, tool, nf_process: nfgen.NFBase):
        inputs = []
        for o in nf_process.outputs:
            # Always use 'val' qualifier
            inp = nfgen.ProcessInput(qualifier=nfgen.InputProcessQualifier.val, name=nf_process.name + o.name)
            inputs.append(inp)

        tool_outputs = cls.prepare_tool_output(tool)

        return inputs, tool_outputs

    @classmethod
    def init_output_process(cls, inputs: List[nfgen.ProcessInput], tool_outputs: Dict):
        script = ""
        for key, val in tool_outputs.items():
            script += f"echo {key}={val} >> {cls.OUTPUT_METADATA_FILENAME}\n"

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
    def init_helper_functions_import(cls, path: Optional[str] = None):
        if path is None:
            path = os.path.join(".", cls.DIR_TOOLS)

        items = [nfgen.ImportItem(name=f.name) for f in cls.generate_generic_functions()]
        return nfgen.Import(items, os.path.join(path, cls.LIB_FILENAME))

    @classmethod
    def init_tool_steps_import(cls, processes: Dict[str, nfgen.Process]) -> List[nfgen.Import]:
        imports = []
        for filename, p in processes.items():
            item = nfgen.ImportItem(name=p.name)
            imp = nfgen.Import([item], os.path.join(".", cls.DIR_TOOLS, filename))
            imports.append(imp)

        return imports

    @classmethod
    def init_process_command_tool(cls, tool: CommandTool, name: Optional[str] = None, provided_inputs: Optional = None) -> nfgen.Process:
        inputs: List[ToolInput] = []
        outputs: List[ToolOutput] = tool.outputs()

        # construct script
        for i in tool.inputs():
            if provided_inputs is not None:
                optional = i.input_type.optional is None or i.input_type.optional is True
                if i.id() in provided_inputs or i.default is not None or not optional:
                    inputs.append(i)
                if isinstance(i.input_type, Filename):
                    inputs.append(i)
            else:
                inputs.append(i)

        script = cls.prepare_script_for_command_tool(tool, inputs)
        pre_script += cls.prepare_expression_inputs(tool, inputs)
        pre_script += cls.prepare_input_vars(tool, inputs)

        resources_var, resource_var_names = cls.prepare_resources_var(tool, name)
        pre_script += resources_var

        pre_script += cls.prepare_inputs_in_selector(tool, inputs, resource_var_names)

        # Logger.debug("INPUT_IN_SELECTORS")
        # Logger.debug(cls.INPUT_IN_SELECTORS)

        process = nfgen.Process(
            name=name or tool.id(),
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            pre_script=pre_script,
        )

        for i in inputs:
            qual = get_input_qualifier_for_inptype(i.input_type)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())

            if isinstance(i.input_type, File):
                # inp.as_process_param = f"Channel.fromPath({cls.PARAM_VAR}).collect()"
                inp.as_process_param = cls.CHANNEL_PARAM

            process.inputs.append(inp)

        for o in outputs:
            output_type = o.output_type
            selector = o.selector

            qual = get_output_qualifier_for_outtype(output_type)
            expression = cls.unwrap_expression(selector, inputs_dict=tool.inputs_map(), tool=tool, for_output=True)

            # #TODO: make this tidier
            if isinstance(output_type, Array):
                if isinstance(output_type.subtype(), (File, Directory)):
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
                is_optional=output_type.optional
            )
            process.outputs.append(out)

        return process

    @classmethod
    def init_process_python_code_tool(cls, tool: CodeTool, name: Optional[str] = None,
                                      provided_inputs: Optional = None) -> nfgen.Process:
        inputs: List[TInput] = []
        outputs: List[TOutput] = tool.outputs()

        # construct script
        for i in tool.inputs():
            if provided_inputs is not None:
                optional = i.intype.optional is None or i.intype.optional is True
                if i.id() in provided_inputs or i.default is not None or not optional:
                    inputs.append(i)
                #TODO: handle input type that works more like arguments
                # note: is Filename set to optional=True, but they have a default value that is not set to i.default
                if isinstance(i.intype, Filename):
                    inputs.append(i)
            else:
                inputs.append(i)

        script = cls.prepare_script_for_python_code_tool(tool, inputs)
        # pre_script += cls.prepare_expression_inputs(tool, inputs)
        # pre_script += cls.prepare_input_vars(tool, inputs)
        #
        # resources_var, resource_var_names = cls.prepare_resources_var(tool, name)
        # pre_script += resources_var
        #
        # pre_script += cls.prepare_inputs_in_selector(tool, inputs, resource_var_names)

        # Logger.debug("INPUT_IN_SELECTORS")
        # Logger.debug(cls.INPUT_IN_SELECTORS)

        process = nfgen.Process(
            name=name or tool.id(),
            script=script,
            script_type=nfgen.ProcessScriptType.script,
        )

        python_file_input = nfgen.ProcessInput(qualifier=nfgen.InputProcessQualifier.path,
                                               name=cls.PYTHON_CODE_FILE_PATH_PARAM.strip("%"),
                                               as_process_param=cls.PYTHON_CODE_FILE_PATH_PARAM)
        process.inputs.append(python_file_input)
        for i in inputs:
            qual = get_input_qualifier_for_inptype(i.intype)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())

            if isinstance(i.intype, File) or \
                    (isinstance(i.intype, Array) and isinstance(i.intype.subtype(), File)):
                # inp.as_process_param = f"Channel.fromPath({cls.PARAM_VAR}).collect()"
                inp.as_process_param = cls.CHANNEL_PARAM

            process.inputs.append(inp)

        for o in outputs:
            qual = nfgen.OutputProcessQualifier.val
            expression = f"file(\"$workDir/{cls.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{o.tag}\").text.replace('[', '').replace(']', '').split(', ')"

            out = nfgen.ProcessOutput(
                qualifier=qual,
                name=o.id(),
                expression=expression
            )
            process.outputs.append(out)


        # for o in outputs:
        #     output_type = o.output_type
        #     selector = o.selector
        #
        #     qual = get_output_qualifier_for_outtype(output_type)
        #     expression = cls.unwrap_expression(selector, inputs_dict=tool.inputs_map(), tool=tool, for_output=True)
        #
        #     # #TODO: make this tidier
        #     if isinstance(output_type, Array):
        #         if isinstance(output_type.subtype(), (File, Directory)):
        #             sub_qual = nfgen.OutputProcessQualifier.path
        #         else:
        #             sub_qual = nfgen.OutputProcessQualifier.val
        #
        #         tuple_elements = expression.strip("][").split(",")
        #         formatted_list = []
        #         for expression in tuple_elements:
        #             sub_exp = nfgen.TupleElementForOutput(qualifier=sub_qual, expression=expression)
        #             formatted_list.append(sub_exp.get_string())
        #
        #         expression = ", ".join(formatted_list)
        #
        #     out = nfgen.ProcessOutput(
        #         qualifier=qual,
        #         name=o.id(),
        #         expression=expression,
        #         is_optional=output_type.optional
        #     )
        #     process.outputs.append(out)

        return process

    @classmethod
    def handle_nf_process_args(cls, nf_process: nfgen.Process, inputsdict: Dict[str, ToolInput], provided_inputs: Dict[str, Any]) -> str:
        args_list = []
        for i in nf_process.inputs:
            if i.name in provided_inputs:
                p = provided_inputs[i.name]
                if p is None:
                    p = "''"
                elif p.startswith("$"):
                    p = p.replace("$", "")
                else:
                    p = f"'{p}'"
            else:
                p = f"'{inputsdict[i.name].default}'" or "''"

            if i.as_process_param:
                # p = i.as_process_param.replace(cls.PARAM_VAR, p)
                # Note: only need to do this for string type input (directly from json file)
                if p.startswith("params."):
                    p = i.as_process_param.replace(cls.CHANNEL_PARAM, f"Channel.fromPath({p}).collect()")

            args_list.append(p)

        args = ", ".join(args_list)

        return args

    @classmethod
    def handle_nf_workflow_args(cls, nf_workflow: nfgen.Workflow, inputsdict: Dict[str, ToolInput], provided_inputs: Dict[str, Any]) -> str:
        pass

    @classmethod
    def init_workflow(cls, workflow, nf_items: Dict[str, Union[nfgen.Process, nfgen.Workflow]], nf_workflow_name: str = ''):
        main = []

        # inputsdict = workflow.inputs_map()
        step_keys = list(workflow.step_nodes.keys())

        Logger.info("nf_processes")
        Logger.info(nf_items.keys())

        for p in nf_items.values():
            Logger.info("process inputs " + p.name)
            for i in p.inputs:
                Logger.info(i.name)

        for step_id in workflow.step_nodes:
            tool = workflow.step_nodes[step_id].tool
            provided_inputs = cls.generate_wf_tool_inputs(tool, step_keys)
            inputsdict = tool.inputs_map()

            # Logger.info(tool.versioned_id())
            # Logger.info("inputsdict")
            # Logger.info(inputsdict)
            # Logger.info(inputsdict.keys())

            #TODO: fetch process or a subworkflow
            item = nf_items.get(tool.versioned_id())

            if isinstance(item, nfgen.Process):
                args = cls.handle_nf_process_args(item, inputsdict, provided_inputs)
            elif isinstance(item, nfgen.Workflow):
                args = cls.handle_nf_workflow_args(item, inputsdict, provided_inputs)

            main.append(f"{item.name}({args})")

        # calling outputs process for Janis to be able to find output files
        args_list = ", ".join([val for val in cls.generate_wf_tool_outputs(workflow).values()])
        main.append(f"outputs({args_list})")

        return nfgen.Workflow(name=nf_workflow_name, main=main)

    @classmethod
    def prepare_output_process_params_for_worfklow(cls, workflow) -> Tuple[List[nfgen.ProcessInput], Dict]:
        wf_outputs = cls.generate_wf_tool_outputs(workflow)

        inputs = []
        outputs = {}
        for key, val in wf_outputs.items():
            inp_var_name = key.replace(".", "")
            # Always use 'val' qualifier
            inp = nfgen.ProcessInput(qualifier=nfgen.InputProcessQualifier.val, name=inp_var_name)
            inputs.append(inp)
            outputs[key] = f"${inp_var_name}"

        return inputs, outputs


    @classmethod
    def init_tool_execution(cls, tool, nf_process: nfgen.Process):
        args_list = []
        for i in nf_process.inputs:

            p = f"params.{i.name}"

            # Extra processing when we need to set up the process input parameters
            if i.as_process_param:
                p = i.as_process_param.replace(cls.CHANNEL_PARAM, f"Channel.fromPath({p}).collect()")

                path_to_python_code_file = posixpath.join("$baseDir", cls.DIR_TOOLS, f"{tool.versioned_id()}.py")
                p = p.replace(cls.PYTHON_CODE_FILE_PATH_PARAM, f"\"{path_to_python_code_file}\"")

            args_list.append(p)

        args = ", ".join(args_list)
        outputs_args = ", ".join(f"{nf_process.name}.out.{o.name}" for o in nf_process.outputs)

        main = [
            f"{nf_process.name}({args})",
            f"outputs({outputs_args})"
        ]

        return nfgen.Workflow(name="", main=main)

    @classmethod
    def generate_wf_tool_inputs(cls, tool: Tool, step_keys: List[str]):
        inputs = {}
        for key in tool.connections:
            if tool.connections[key] is None:
                val = None
            else:
                val = str(tool.connections[key])

                if "inputs." in val:
                    # e.g. replace inputs.fastq to params.
                    val = val.replace("inputs.", "$params.")

                # e.g. replace bwamem.var to bwamem.out.var
                # bwamem.out is nextflow variable to fetch all output from bwamem process
                for tool_id in step_keys:
                    keyword = f"{tool_id}."
                    if keyword in val:
                        val = val.replace(keyword, f"${tool_id}.out.")

            inputs[key] = val

        return inputs

    @classmethod
    def generate_wf_tool_outputs(cls, wf: WorkflowBase):
        step_keys = wf.step_nodes.keys()

        outputs = {}
        for o in wf.output_nodes:
            val = str(wf.output_nodes[o].source)

            if "inputs." in val:
                # e.g. replace inputs.fastq to $fastq
                val = val.replace("inputs.", "params.")

            # e.g. replace bwamem.var to bwamem.out.var
            # bwamem.out is nextflow variable to fetch all output from bwamem process
            for tool_id in step_keys:
                keyword = f"{tool_id}."
                if keyword in val:
                    val = val.replace(keyword, f"{tool_id}.out.")

            outputs[o] = val

        return outputs

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
                    for_output=for_output,
                    in_shell_script=in_shell_script
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
                in_shell_script=in_shell_script,
                tool=tool
            )
        elif isinstance(value, WildcardSelector):
            # raise Exception(
            #     f"A wildcard selector cannot be used as an argument value for '{debugkwargs}' {tool.id()}"
            # )
            return f"'{value.wildcard}'"
        elif isinstance(value, Operator):
            unwrap_expression_wrap = lambda exp: cls.unwrap_expression(
                exp,
                quote_string=quote_string,
                tool=tool,
                for_output=for_output,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                in_shell_script=in_shell_script,
                **debugkwargs,
            )

            return value.to_nextflow(unwrap_expression_wrap, *value.args)

        elif callable(getattr(value, "nextflow", None)):
            return value.Nextflow()

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
                elif hasattr(tinp, "localise_file") and tinp.localise_file:
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
                                 tool,
                                 skip_inputs_lookup=False,
                                 in_shell_script=False,
                                 for_output=False
                                 ):
        # TODO: Consider grabbing "path" of File

        if tool.id() not in cls.INPUT_IN_SELECTORS:
            cls.INPUT_IN_SELECTORS[tool.versioned_id()] = set()

        cls.INPUT_IN_SELECTORS[tool.versioned_id()].add(selector.input_to_select)

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
        # Logger.debug("values_provided_from_tool")
        # Logger.debug(values_provided_from_tool)
        for i in tool.tool_inputs():
            val = ad.get(i.id(), values_provided_from_tool.get(i.id()))

            inputsdict = tool.inputs_map()

            if isinstance(i.intype, Filename):
                val = cls.unwrap_expression(i.intype.generated_filename(), inputs_dict=inputsdict, tool=tool)
            elif isinstance(i.intype, Boolean):
                val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""
                if val == "True":
                    val = True
                if val == "False":
                    val = False
            if isinstance(i.intype, File):
                if hasattr(i.intype, 'secondary_files') and callable(i.intype.secondary_files):
                    if i.intype.secondary_files() is not None:
                        primary_file = val
                        secondary_files = []
                        for suffix in i.intype.secondary_files():
                            sec_file = apply_secondary_file_format_to_filename(primary_file, suffix)
                            secondary_files.append(sec_file)

                        val = [primary_file] + secondary_files

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
    def has_secondary_files(cls, i: ToolInput):
        if hasattr(i.input_type, 'secondary_files') and callable(i.input_type.secondary_files):
            if i.input_type.secondary_files() is not None:
                return True
        return False

    @classmethod
    def prepare_script_for_python_code_tool(cls, tool: CodeTool, inputs):
        PYTHON_SHEBANG = "#!/usr/bin/env python"
        # code_block = tool.prepared_script(SupportedTranslation.Nextflow)
        # code_block = code_block.replace('\\', '\\\\')
        # code_block = code_block.replace('$', '\\$')
        # code_block = code_block.replace('"""', '\\"\\"\\"')

        python_script_filename = f"{tool.versioned_id()}"
        # python_script_path = os.path.join("$baseDir", cls.DIR_TOOLS, python_script_filename)
        # python_script_path = os.path.join("$baseDir", cls.DIR_TOOLS, python_script_filename)


        # TODO: handle args of type list of string (need to quote them)
        # args = ", ".join(i.tag for i in inputs)
        all_args = []
        for i in inputs:
            arg_value = f"${i.tag}"
            if isinstance(i.intype, Array) and not isinstance(i.intype.subtype(), (Int, Float, Double)):
                # arg_value = f"\"{arg_value}\".strip(\"][\").split(\", \")"
                arg_value = f"\"{arg_value}\".split(\" \")"

            elif not isinstance(i.intype, (Array, Int, Float, Double)):
                arg_value = f"\"{arg_value}\""

            all_args.append(arg_value)

        args = ", ".join(a for a in all_args)

        script = f"""{PYTHON_SHEBANG}
from {python_script_filename} import code_block
import os
import json

result = code_block({args})

work_dir = os.getenv("PYENV_DIR")
for key in result:
    with open(os.path.join("$workDir", f"{cls.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{{key}}"), "w") as f:
        f.write(json.dumps(result[key]))
"""

        return script

    @classmethod
    def prepare_script_for_command_tool(cls, tool: CommandTool, inputs):
        bc = tool.base_command()
        pargs = []

        if bc:
            pargs.append(" ".join(bc) if isinstance(bc, list) else str(bc))

        arguments = []
        if isinstance(tool, CommandTool):
            arguments = tool.arguments() or []

        args = [a for a in arguments if a.position is not None or a.prefix is not None]
        args += [a for a in inputs if a.position is not None or a.prefix is not None]

        args = sorted(args, key=lambda a: (a.prefix is None))
        args = sorted(args, key=lambda a: (a.position or 0))

        prefix = "  "
        for a in args:
            if isinstance(a, ToolInput):
                pargs.append(f"${a.id()}WithPrefix")
            elif isinstance(a, ToolArgument):
                inputsdict = tool.inputs_map()
                expression = cls.unwrap_expression(
                    a.value,
                    inputsdict=inputsdict,
                    tool=tool,
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
{main_script} > {nfgen.Process.TOOL_STDOUT_FILENAME}

"""

    @classmethod
    def prepare_tool_output(cls, tool: Tool):
        inputsdict = tool.inputs_map()

        outputs = {}
        for out in tool.outputs():
            if isinstance(out, TOutput):
                output_type = out.outtype
            elif isinstance(out, ToolOutput):
                output_type = out.output_type
            else:
                raise Exception("unknown output object type")


            if isinstance(output_type, Stdout):
                val = "STDOUT"
            elif isinstance(output_type, Stderr):
                val = "STDERR"
            else:
                val = f"${tool.id()}{out.tag}"

            outputs[out.tag] = val

        return outputs

    @classmethod
    def prepare_expression_inputs(cls, tool, inputs):

        inputsdict = tool.inputs_map()

        script_lines = []
        for i in inputs:
            if hasattr(i, "input_type"):
                input_type = i.input_type
            elif hasattr(i, "intype"):
                input_type = i.intype
            else:
                raise Exception("Failed to get input type attribute")

            if isinstance(input_type, Filename):
                val = cls.unwrap_expression(input_type.generated_filename(), inputs_dict=inputsdict, tool=tool)

                code = f"""
def {i.id()} = {val}
"""
                script_lines.append(code)

        return "".join(script_lines)

    @classmethod
    def prepare_input_vars(self, tool, inputs):
        pre_script_lines = []
        for a in inputs:
            arg_name = ""
            if isinstance(a, ToolInput):
                arg_name = a.id()
            elif isinstance(a, ToolArgument):
                continue
            else:
                raise Exception("unknown input type")

            if isinstance(a.input_type, Array):
                arg_value = f"{arg_name}.join(' ')"
            elif isinstance(a.input_type, (File)) and self.has_secondary_files(a):
                arg_value = f"{arg_name}[0]"
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
    def prepare_resources_var(cls, tool, name: Optional[str] = None):
        pre_script_lines = []
        var_names = []
        for k, v in cls.build_resources_input(
                tool,
                hints=None,
        ).items():
            prefix = ''
            if name is not None:
                prefix = f"{name}_"

            code = f"""
def {k} =  params.{prefix}{k}
"""
            var_names.append(k)
            pre_script_lines.append(code)

        return "".join(pre_script_lines), var_names

    @classmethod
    def prepare_inputs_in_selector(cls, tool, inputs: List[ToolInput], resources: List[str]):
        pre_script_lines = []

        if tool.versioned_id() not in cls.INPUT_IN_SELECTORS:
            return ""

        input_keys = [i.id() for i in inputs]

        for k in cls.INPUT_IN_SELECTORS[tool.versioned_id()]:
            if k not in input_keys and k not in resources:
                val = "''"

                code = f"""
def {k} =  {val}
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
