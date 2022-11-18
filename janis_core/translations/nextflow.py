from copy import deepcopy
import os
import re
import posixpath
from typing import Tuple, Dict, List, Optional, Union, Any
from janis_core.translations.nfgen import format_process_call
from janis_core.types import (
    DataType,
    Array,
    File,
    Int,
    Float,
    Double,
    Directory,
    Stdout,
    Filename,
)

from janis_core.utils.scatter import ScatterDescription, ScatterMethod

from janis_core.tool.commandtool import (
    CommandTool,
    ToolInput,
    ToolOutput,
    Tool,
    TOutput,
    TInput,
)
from janis_core.code.codetool import CodeTool
from janis_core.code.pythontool import PythonTool
from janis_core.translations.translationbase import TranslatorBase
from janis_core.workflow.workflow import StepNode, Workflow, WorkflowBase
from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.translations import nfgen
from janis_core.translations.nfgen import settings

# class methods dont make sense. dislike this approach. 
# whole approach is far too complex. no need for oop. 
# should have just defined a translation interface (protocol class in python) 
# specifying the method signatures each translator needs.

# class File:
#     ...

# def translate_workflow(workflow: Workflow) -> list[File]:
#     """translates a janis workflow to nextflow format."""
#     ...

# def translate_tool(tool: Tool) -> list[File]:
#     """translates a janis tool to nextflow"""
#     ...


class NextflowTranslator(TranslatorBase):
    INPUT_IN_SELECTORS: dict[str, Any] = {}

    def __init__(self):
        super().__init__(name="nextflow")

    @classmethod
    def translate_workflow(
        cls,
        jworkflow: Workflow,
        with_container: bool = True,  
        with_resource_overrides: bool = False,
        allow_empty_container: bool = False,
        container_override: Optional[dict[str, str]] = None,
        render_comments: bool = True
    ) -> Tuple[Any, dict[str, Any]]:

        scope: list[str] = []
        files: dict[str, nfgen.NFFile] = {}
        processes: dict[str, nfgen.Process] = {}
        subworkflows: dict[str, nfgen.Workflow] = {}
        operations: dict[str, nfgen.ChannelOperation] = {}

        # register params and channels for workflow inputs
        nfgen.channels.register(workflow=jworkflow)
        nfgen.params.register(the_entity=jworkflow, scope=scope)

        # parse each step to a NFFile
        for step in jworkflow.step_nodes.values():
            current_scope = deepcopy(scope)
            current_scope.append(step.identifier)         
            nf_items = cls.gen_items_for_step(
                step,
                scope=current_scope,
                with_container=with_container,
                with_resource_overrides=with_resource_overrides,
                allow_empty_container=allow_empty_container,
                container_override=container_override,
            )
            
            nf_file = nfgen.NFFile(
                #imports=[cls.init_helper_functions_import(".")], 
                imports=[], 
                items=nf_items,
                #name=f"{step.id()}_{step.tool.versioned_id()}"
                name=step.id()
            )

            files[nf_file.name] = nf_file
            for item in nf_items:
                if isinstance(item, nfgen.Process):
                    processes[step.id()] = item
                elif isinstance(item, nfgen.Workflow):
                    subworkflows[step.id()] = item
                elif isinstance(item, nfgen.ChannelOperation):
                    operations[step.id()] = item
                else:
                    raise NotImplementedError

        # handle imports
        imports = cls.collect_workflow_imports(files)
        
        # create outputs collection process
        # out_process_inp, out_process_out = cls.prepare_output_process_params_for_workflow(jworkflow)
        # output_p = cls.gen_output_process(out_process_inp, out_process_out)
        
        # create object & NFFile for workflow
        workflow = cls.gen_workflow(
            janis=jworkflow, 
            processes=processes,
            operations=operations,
            subworkflows=subworkflows
        )
        channels = nfgen.channels.channel_register
        nf_file = nfgen.NFFile(imports=imports, items=[channels, workflow])

        # generate strings for each file
        tool_scripts: dict[str, str] = {name: nffile.get_string() for name, nffile in files.items()}
        return (nf_file.get_string(), tool_scripts)
        
    @classmethod
    def gen_items_for_step(
        cls,
        step: StepNode,
        scope: list[str],
        # step_keys: Optional[list[str]] = None,
        with_container: bool = True, 
        with_resource_overrides: bool = False,
        allow_empty_container: bool = False,
        container_override: Optional[dict[str, str]] = None,
    ) -> list[nfgen.Process | nfgen.Workflow | nfgen.ChannelOperation]:
        """
        For each of the workflow step, we need to generate a Nextflow subworkflow or process object

        :param tool:
        :type tool:
        :param step_id:
        :type step_id:
        :param step_keys:
        :type step_keys:
        :param with_container:
        :type with_container:
        :param with_resource_overrides:
        :type with_resource_overrides:
        :param allow_empty_container:
        :type allow_empty_container:
        :param container_override:
        :type container_override:
        :return:
        :rtype:
        """

        nf_items: list[nfgen.Process | nfgen.Workflow | nfgen.ChannelOperation] = []
        if isinstance(step.tool, CommandTool):
            if step.scatter and step.scatter.method == ScatterMethod.cross:
                operation = nfgen.channels.gen_scatter_cross_operation(step.sources, step.scatter)
                nf_items.append(operation)
            process = cls.gen_process_from_cmdtool(step.tool, step.sources, scope)
            process = cls.handle_container(
                step.tool, process, with_container, allow_empty_container, container_override
            )
            nf_items.append(process)
            return nf_items

        elif isinstance(step.tool, PythonTool):
            process = cls.gen_process_from_codetool(step.tool, step.sources, scope)
            process = cls.handle_container(
                step.tool, process, with_container, allow_empty_container, container_override
            )
            return [process]

        elif isinstance(step.tool, WorkflowBase):
            subworkflow = step.tool
            # register params for subworkflow inputs.
            # no channels to register for subworkflow.
            nfgen.params.register(the_entity=subworkflow, sources=step.sources, scope=scope)

            for substep in subworkflow.step_nodes.values():
                current_scope = deepcopy(scope)
                current_scope.append(step.identifier)
                nf_items += cls.gen_items_for_step(
                    substep,
                    current_scope,
                    # step_keys=sub_step_keys,
                    with_container=with_container,
                    with_resource_overrides=with_resource_overrides,
                    allow_empty_container=allow_empty_container,
                    container_override=container_override,
                )
            nf_items = [cls.gen_subworkflow(step, step.id(), nf_items)] + nf_items
            return nf_items

        elif isinstance(step.tool, CodeTool):
            raise Exception("Only PythonTool code tool is supported for the moment")
        
        else:
            raise Exception(
                f"Nextflow translation for this {step.tool.__class__} is not yet supported"
            )

    @classmethod
    def translate_tool_internal(
        cls,
        tool: CommandTool,
        with_container: bool = True,
        with_resource_overrides: bool = False,
        allow_empty_container: bool = False,
        container_override: Optional[dict[str, str]] = None,
        render_comments: bool = True
    ) -> str:
        """
        Generate Nextflow translation for Janis command line tool
        main_p - the main tool process 
        output_p - the process to collect outputs

        :param tool:
        :type tool:
        :param with_container:
        :type with_container:
        :param with_resource_overrides:
        :type with_resource_overrides:
        :param allow_empty_container:
        :type allow_empty_container:
        :param container_override:
        :type container_override:
        :return:
        :rtype:
        """
        settings.MODE = 'tool'
        file_items: list[nfgen.NFBase] = []
        scope: list[str] = []
        values: dict[str, Any] = {}

        main_p = cls.gen_process_from_cmdtool(tool, values, scope)
        main_p = cls.handle_container(
            tool, main_p, with_container, allow_empty_container, container_override
        )
        # output_p_inp, output_p_out = cls.prepare_output_process_params_for_tool(tool, main_p)
        # output_p = cls.gen_output_process(inputs=output_p_inp, tool_outputs=output_p_out)
        workflow = cls.gen_process_workflow(tool, process=main_p)

        file_items.append(param_block)
        #file_items.append(channel_block)
        file_items.append(main_p)
        #file_items.append(output_p)
        file_items.append(workflow)

        nf_file = nfgen.NFFile(
            #imports=[cls.init_helper_functions_import('.')], 
            imports=[], 
            items=file_items
        )
        return nf_file.get_string()

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool: CodeTool,
        with_container: bool = True,  
        allow_empty_container: bool = False,
        container_override: Optional[dict[str, str]] = None,
        render_comments: bool = True
    ) -> str:
        """
        Generate Nextflow translation for Janis code tool

        :param tool:
        :type tool:
        :param with_container:
        :type with_container:
        :param allow_empty_container:
        :type allow_empty_container:
        :param container_override:
        :type container_override:
        :return:
        :rtype:
        """
        raise NotImplementedError
        if isinstance(tool, PythonTool):
            process = cls.gen_process_from_codetool(tool=tool, values={}, scope=[])
            process = cls.handle_container(
                tool, process, with_container, allow_empty_container, container_override
            )

            #imports = [cls.init_helper_functions_import()]
            imports = []

            # (
            #     out_process_inp,
            #     out_process_out,
            # ) = cls.prepare_output_process_params_for_tool(tool, process)

            items = [
                process,
                # cls.gen_output_process(
                #     inputs=out_process_inp, tool_outputs=out_process_out
                # ),
                cls.gen_process_workflow(tool, process=process),
            ]
            nf_file = nfgen.NFFile(imports=imports, items=items)

            return nf_file.get_string()
        else:
            raise Exception("Only PythonTool code tool is supported for the moment.")

    @classmethod
    def gen_subworkflow(
        cls, step: StepNode, name: str, nf_items: List[nfgen.NFBase]
    ) -> nfgen.Workflow:
        """
        Generate a Nextflow Workflow object for a Janis subworkflow inside a workflow

        :param subworkflow:
        :type subworkflow:
        :param name:
        :type name:
        :param nf_items:
        :type nf_items:
        :return:
        :rtype:
        """
        body: list[str] = []
        take = []
        emit = []
        subworkflow: Workflow = step.tool

        inputsdict = subworkflow.inputs_map()
        step_keys = list(subworkflow.step_nodes.keys())
        step_inputs = cls.gen_step_inval_dict(step.tool, step.sources)

        for i in subworkflow.input_nodes:
            if i not in step_inputs:
                if subworkflow.input_nodes[i].default is not None:
                    step_inputs[i] = subworkflow.input_nodes[i].default

        for key in subworkflow.connections:
            as_param = None
            input_type = inputsdict.get(key).intype
            if isinstance(input_type, File):
                # inp.as_process_param = f"Channel.fromPath({settings.PARAM_VAR}).collect()"
                as_param = settings.LIST_OF_FILES_PARAM
            elif isinstance(input_type, Array) and isinstance(
                input_type.subtype(), Array
            ):
                as_param = settings.LIST_OF_FILE_PAIRS_PARAM

            take.append(nfgen.WorkflowInput(name=key, as_param=as_param))

        for step_id in subworkflow.step_nodes:
            nf_item = [i for i in nf_items if i.name == f"{name}_{step_id}"][0]
            tool = subworkflow.step_nodes[step_id].tool

            step_inputs = cls.gen_step_inval_dict(
                step.tool,
                step.sources,
                inputs_replacement="$",
                tool_id_prefix=f"${name}_",
            )

            step_inputs = cls.apply_outer_workflow_inputs(
                subworkflow, tool, step_inputs, step_inputs
            )

            args = cls.handle_process_args(
                tool,
                nf_item.inputs,
                step_inputs,
                input_param_prefix=f"{name}_{step_id}_",
                step_key_prefix=f"{name}_",
                workflow=subworkflow,
                scatter=subworkflow.step_nodes[step_id].scatter,
            )

            body.append(format_process_call(f'{name}_{step_id}', args))

        wf_outputs = cls.gen_wf_tool_outputs(subworkflow, f"{name}_")
        for o in wf_outputs:
            emit.append(nfgen.WorkflowOutput(name=o, expression=wf_outputs[o]))

        return nfgen.Workflow(name=name, main=body, take=take, emit=emit)

    @classmethod
    def gen_workflow(
        cls,
        janis: Workflow,
        processes: dict[str, nfgen.Process],
        operations: dict[str, nfgen.ChannelOperation],
        subworkflows: dict[str, nfgen.Workflow],
        nf_workflow_name: str = "",
    ) -> nfgen.Workflow:
        """
        Generate a Nextflow Workflow object

        :param workflow:
        :type workflow:
        :param nf_items:
        :type nf_items:
        :param nf_workflow_name:
        :type nf_workflow_name:
        :return:
        :rtype:
        """
        body: list[str] = []

        for step_id, step in janis.step_nodes.items():
            tool = janis.step_nodes[step_id].tool
            step_inputs = cls.gen_step_inval_dict(step.tool, step.sources)

            # if there are operations, add these to body before the process/subworkflow call
            if step_id in operations:
                operation = operations[step_id]
                body.append(operation.get_string())

            if step_id in processes:
                entity_name = processes[step_id].name
                entity_inputs = processes[step_id].inputs
            elif step_id in subworkflows:            
                entity_name = subworkflows[step_id].name
                entity_inputs = subworkflows[step_id].take
            else:
                raise NotImplementedError

            args = cls.handle_process_args(
                tool,
                entity_inputs,
                step_inputs,
                input_param_prefix=f"{step_id}_",
                workflow=janis,
                scatter=janis.step_nodes[step_id].scatter,
            )
            body.append(format_process_call(entity_name, args))

        # calling outputs process for Janis to be able to find output files
        # args_list = [val for val in cls.gen_wf_tool_outputs(janis).values()]
        # body.append(format_process_call(settings.FINAL_STEP_NAME, args_list))
        return nfgen.Workflow(name=nf_workflow_name, main=body)

    @classmethod
    def apply_outer_workflow_inputs(
        cls,
        workflow: WorkflowBase,
        tool: Tool,
        provided_inputs: Dict[str, Any],
        wf_provided_inputs: Dict[str, Any],
    ) -> dict:
        """
        Generate a dictionary of user provided inputs for a tool by
        taking into considerations inputs provided to its parent workflow

        :param workflow:
        :type workflow:
        :param provided_inputs: inputs provided to this tool
        :type provided_inputs:
        :param wf_provided_inputs: inputs provided to the outer workflow this tool is called from
        :type wf_provided_inputs:
        :return:
        :rtype:
        """
        # Apply value from workflow that is calling it
        wf_input_names = workflow.connections.keys()
        step_keys = workflow.step_nodes.keys()
        wf_inputsdict = workflow.inputs_map()
        tool_inputsdict = tool.inputs_map()

        for key in provided_inputs:
            val = provided_inputs[key]

            # Here, we are looking for input variables that are not an input of the subworkflow
            # but, it is pointing to a variable of the outer workflow that calls this workflow
            if (
                key not in wf_input_names
                and isinstance(val, str)
                and val.startswith("$")
            ):
                if key in wf_provided_inputs:
                    provided_inputs[key] = str(wf_provided_inputs[key])

            # if we are just using the variable name from the nextflow subworkflow
            for wf_input in wf_inputsdict:
                if val is not None and val == f"${wf_input}":
                    wf_inp_type = wf_inputsdict.get(wf_input).intype
                    tool_inp_type = (
                        tool_inputsdict.get(key).intype
                        or tool_inputsdict.get(key).input_type
                    )
                    if (
                        isinstance(wf_inp_type, File)
                        and wf_inp_type.has_secondary_files()
                    ):
                        if (
                            tool_inp_type.is_base_type(File)
                            and not tool_inp_type.has_secondary_files()
                        ) or (
                            isinstance(tool_inp_type, Array)
                            and isinstance(tool_inp_type.subtype(), File)
                            and not tool_inp_type.subtype().has_secondary_files()
                        ):
                            provided_inputs[key] += ".map{ tuple -> tuple[0] }"

            for step_key in step_keys:
                if isinstance(val, str) and f"{step_key}.out." in val:

                    parts = val.strip("][").split(f"{step_key}.out.")
                    step_output_var = parts[1]
                    step_outputs = workflow.step_nodes.get(step_key).outputs()

                    src_type = step_outputs.get(step_output_var).outtype
                    tool_inp_type = (
                        tool_inputsdict.get(key).intype
                        or tool_inputsdict.get(key).input_type
                    )

                    # any file types or array of files
                    if src_type.is_base_type(File) and src_type.has_secondary_files():
                        if (
                            tool_inp_type.is_base_type(File)
                            and not tool_inp_type.has_secondary_files()
                        ) or (
                            isinstance(tool_inp_type, Array)
                            and isinstance(tool_inp_type.subtype(), File)
                            and not tool_inp_type.subtype().has_secondary_files()
                        ):
                            provided_inputs[key] += ".map{ tuple -> tuple[0] }"

        # Now, get rid of anything that points to workflow input variables but is not provided as workflow inputs
        for key, val in provided_inputs.items():
            if isinstance(val, str) and val.startswith("$"):
                val = val.strip("$")
                if val in wf_inputsdict and val not in wf_provided_inputs:
                    provided_inputs[key] = None

        return provided_inputs

    @classmethod
    def handle_container(
        cls,
        tool: Tool,
        process: nfgen.Process,
        with_container: bool = True,
        allow_empty_container: bool = False,
        container_override: Optional[dict[str, str]] = None,
    ) -> nfgen.Process:
        """
        Add container information to a Nexflow Process object

        :param tool:
        :type tool:
        :param process:
        :type process:
        :param with_container:
        :type with_container:
        :param allow_empty_container:
        :type allow_empty_container:
        :param container_override:
        :type container_override:
        :return:
        :rtype:
        """
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
                    container=cls.unwrap_expression(container, quote_string=False, tool=tool)
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
    def gen_python_script(cls, tool: PythonTool) -> str:
        """
        Generate python script to be included inside a Nextflow process

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        return tool.prepared_script(SupportedTranslation.NEXTFLOW)

    @classmethod
    def translate_helper_files(cls, tool) -> Dict[str, str]:
        """
        Generate a dictionary of helper files to run Nextflow.
        Key of the dictionary is the filename, the value is the file content

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        helpers = {}

        # Groovy helper functions
        lib_file = nfgen.NFFile(imports=[], items=cls.gen_generic_functions())
        helpers[settings.LIB_FILENAME] = lib_file.get_string()

        helpers = cls.gen_python_code_files(tool, helpers)

        return helpers

    @classmethod
    def gen_python_code_files(cls, tool: PythonTool, helpers: dict):
        # Python files for Python code tools
        if isinstance(tool, PythonTool):
            helpers["__init__.py"] = ""
            #helpers[f"{tool.versioned_id()}.py"] = cls.gen_python_script(tool)
            helpers[f"{tool.id()}.py"] = cls.gen_python_script(tool)
            return helpers

        elif isinstance(tool, WorkflowBase):
            for step_id in tool.step_nodes:
                step_tool = tool.step_nodes[step_id].tool
                helpers = cls.gen_python_code_files(step_tool, helpers)

        return helpers

    # @classmethod
    # def prepare_output_process_params_for_tool(
    #     cls, tool: Tool, nf_process: nfgen.NFBase
    # ):
    #     """
    #     Every one of our tools will call a Nextflow process named "janis_outputs".
    #     This process will collect all the final outputs for this tool or workflow.
    #     This allows Janis to process outputs more easily.

    #     This function is used to generate the inputs and outputs for "janis_outputs" Nextflow process.

    #     :param tool:
    #     :type tool:
    #     :param nf_process: This is the Nextflow Process of the actual tool
    #     :type nf_process:
    #     :return:
    #     :rtype:
    #     """
    #     inputs = []
    #     for o in nf_process.outputs:
    #         # Always use 'val' qualifier
    #         inp = nfgen.ProcessInput(
    #             qualifier=nfgen.InputProcessQualifier.val, name=nf_process.name + o.name
    #         )
    #         inputs.append(inp)

    #     tool_outputs = cls.prepare_tool_output(tool)

    #     return inputs, tool_outputs

    # @classmethod
    # def prepare_output_process_params_for_workflow(
    #     cls, workflow: WorkflowBase
    # ) -> Tuple[List[nfgen.ProcessInput], dict[str, Any]]:
    #     """
    #     Every one of our tools will call a Nextflow process named "janis_outputs".
    #     This process will collect all the final outputs for this tool or workflow.
    #     This allows Janis to process outputs more easily.

    #     This function is used to generate the inputs and outputs for "janis_outputs" Nextflow workflow.

    #     :param workflow:
    #     :type workflow:
    #     :return:
    #     :rtype:
    #     """
    #     wf_outputs = cls.gen_wf_tool_outputs(workflow)
    #     output_dict = workflow.outputs_map()

    #     inputs = []
    #     outputs = {}

    #     for key, val in wf_outputs.items():
    #         inp_var_name = key.replace(".", "")
    #         # Always use 'val' qualifier
    #         inp = nfgen.ProcessInput(
    #             qualifier=nfgen.InputProcessQualifier.val, name=inp_var_name
    #         )
    #         inputs.append(inp)

    #         output_var = inp_var_name
    #         if key in output_dict:
    #             if (
    #                 output_dict[key].outtype.is_array()
    #                 and isinstance(output_dict[key].outtype.subtype(), File)
    #                 and output_dict[key].outtype.subtype().has_secondary_files()
    #             ):
    #                 output_var = f"{inp_var_name}.map{{ item -> item[0] }}"
    #             elif (
    #                 isinstance(output_dict[key].outtype, File)
    #                 and output_dict[key].outtype.has_secondary_files()
    #             ):
    #                 output_var = f"{inp_var_name}[0]"

    #         outputs[key] = f"${{{output_var}}}"

    #     return inputs, outputs

    # @classmethod
    # def gen_output_process(
    #     cls, inputs: List[nfgen.ProcessInput], tool_outputs: dict[str, Any]
    # ) -> nfgen.Process:
    #     """
    #     Every one of our tools will call a Nextflow process named "janis_outputs".
    #     This process will collect all the final outputs for this tool or workflow.
    #     This allows Janis to process outputs more easily.

    #     This function creates this "janis_outputs" Nextflow process.

    #     :param inputs:
    #     :type inputs:
    #     :param tool_outputs:
    #     :type tool_outputs:
    #     :return:
    #     :rtype:
    #     """

    #     # The script is simply adding tool outputs as key-value pairs into a text file
    #     script = ""
    #     for key, val in tool_outputs.items():
    #         script += f"echo {key}={val} >> {settings.OUTPUT_METADATA_FILENAME}\n"

    #     outputs = [
    #         nfgen.ProcessOutput(
    #             qualifier=nfgen.OutputProcessQualifier.path,
    #             name="janis_output_metadata",
    #             expression=f"'{settings.OUTPUT_METADATA_FILENAME}'",
    #         )
    #     ]

    #     process = nfgen.Process(
    #         name=settings.FINAL_STEP_NAME,
    #         script=script,
    #         script_type=nfgen.ProcessScriptType.script,
    #         inputs=inputs,
    #         outputs=outputs,
    #     )

    #     process.directives.append(nfgen.CacheDirective(enabled=False))

    #     return process

    @classmethod
    def init_helper_functions_import(cls, path: Optional[str] = None) -> nfgen.Import:
        """
        Generate Nextflow Import object to import our Groovy helper functions

        :param path:
        :type path:
        :return:
        :rtype:
        """
        if path is None:
            path = os.path.join(".", cls.DIR_TOOLS)


        items: list[nfgen.ImportItem] = []

        # TODO logic here to decide which helper function imports are needed
        # if isinstance(inp.input_type, Array):

        items = [
            nfgen.ImportItem(name=f.name) for f in cls.gen_generic_functions()
        ]

        

        return nfgen.Import(items, os.path.join(path, settings.LIB_FILENAME))

    @classmethod
    def collect_workflow_imports(
        cls, files: dict[str, nfgen.NFFile]
    ) -> list[nfgen.Import]:
        """
        Generate a list of Nextflow Import objects for each of the steps in a workflow.

        :param files:
        :type files:
        :return:
        :rtype:
        """
        imports: list[nfgen.Import] = []
        #imports += [cls.init_helper_functions_import()]
        
        for filename, nf_file in files.items():
            file_imports = []
            for process in nf_file.items:
                i_item = nfgen.ImportItem(name=process.name)
                file_imports.append(i_item)

            imp = nfgen.Import(
                file_imports, os.path.join(".", cls.DIR_TOOLS, filename)
            )
            imports.append(imp)

        return imports

    @classmethod
    def gen_process_from_cmdtool(
        cls,
        tool: CommandTool,
        values: dict[str, Any],   # values fed to tool inputs (step translation)
        scope: list[str],
    ) -> nfgen.Process:
        """
        Generate a Nextflow Process object for a Janis Command line tool

        :param tool:
        :type tool:
        :param name: Generally, this is a workflow step id, so that we can prefix variables or process names
        :type name:
        :param provided_inputs:
        :type provided_inputs:
        :return:
        :rtype:
        """

        # nfgen.params.register(the_entity=tool, sources=values, scope=scope)

        process_name = scope[-1] if scope else tool.id()
        pre_script, script = nfgen.gen_script_for_cmdtool(
            tool=tool,
            scope=scope,
            values=values,
            input_in_selectors=cls.INPUT_IN_SELECTORS, 
            stdout_filename=settings.TOOL_STDOUT_FILENAME
        )
 
        process = nfgen.Process(
            name=process_name,
            pre_script=pre_script,
            script=script,
            script_type=nfgen.ProcessScriptType.script,
        )
        
        process_ids = nfgen.utils.get_process_input_ids(tool, values)
        process_inputs = nfgen.utils.items_with_id(tool.inputs(), process_ids)
        for i in process_inputs:
            qual = cls.get_input_qualifier_for_inptype(i.input_type)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())

            process.inputs.append(inp)

        #resources = cls.build_resources_input(tool, hints=None)
        resources = {}
        process.outputs = cls.gen_outputs_for_process(process_name, tool)  # TODO change to scope
        process.directives = cls.gen_directives_for_process(
            resources, scope
        )

        return process

    @classmethod
    def gen_directives_for_process(
        cls, resources: dict[str, Any], scope: Optional[list[str]]=None
    ) -> list[nfgen.ProcessDirective]:
        scope = scope if scope else []

        nf_directives: list[nfgen.ProcessDirective] = []
        nf_directives.append(nfgen.PublishDirDirective(scope))

        # Add directives for input resources
        for res, val in resources.items():
            if res.endswith("runtime_cpu"):
                nf_directives.append(nfgen.CpusDirective(scope, res, val))
            elif res.endswith("runtime_memory"):
                nf_directives.append(nfgen.MemoryDirective(scope, res, val))
            elif res.endswith("runtime_seconds"):
                nf_directives.append(nfgen.TimeDirective(scope, res, val))
            elif res.endswith("runtime_disk"):
                nf_directives.append(nfgen.DiskDirective(scope, res, val))
        return nf_directives

    @classmethod
    def gen_outputs_for_process(
        cls, process_name: str, tool: Tool
    ) -> List[nfgen.ProcessOutput]:
        """
        Generate a list of tool outputs in the form of nfgen.ProcessOutput objects

        :param process_name:
        :type process_name:
        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        outputs: List[ToolOutput] = tool.outputs()

        nf_outputs = []
        for o in outputs:
            output_type = o.output_type
            selector = o.selector

            qual = cls.get_output_qualifier_for_outtype(output_type)
            expression = cls.unwrap_expression(
                selector, inputs_dict=tool.inputs_map(), tool=tool, for_output=True
            )

            if isinstance(output_type, Array):
                if isinstance(output_type.subtype(), (File, Directory)):
                    sub_qual = nfgen.OutputProcessQualifier.path
                else:
                    sub_qual = nfgen.OutputProcessQualifier.val

                tuple_elements = expression.strip("][").split(",")
                formatted_list = []
                for expression in tuple_elements:
                    sub_exp = nfgen.TupleElementForOutput(
                        qualifier=sub_qual, expression=expression
                    )
                    formatted_list.append(sub_exp.get_string())

                expression = ", ".join(formatted_list)
            elif isinstance(output_type, Stdout):
                expression = f'"{settings.TOOL_STDOUT_FILENAME}_{process_name}"'
            elif isinstance(output_type, File) and output_type.has_secondary_files():
                sub_qual = nfgen.OutputProcessQualifier.path
                tuple_elements = [expression]

                primary_ext = output_type.extension
                secondary_ext = []

                if o.secondaries_present_as is not None:
                    secondaries_present_as = o.secondaries_present_as
                else:
                    secondaries_present_as = {}

                for ext in output_type.secondary_files():
                    if ext in secondaries_present_as:
                        secondary_ext.append(secondaries_present_as[ext])
                    else:
                        secondary_ext.append(ext)

                for ext in secondary_ext:
                    replacement = primary_ext + ext
                    if ext.startswith("^"):
                        replacement = ext[1:]

                    sec_exp = None
                    if primary_ext in expression:
                        sec_exp = expression.replace(primary_ext, replacement)
                    elif ".name" in expression:
                        sec_exp = expression.replace(
                            ".name", f".baseName + '{replacement}'"
                        )

                    if sec_exp is not None:
                        tuple_elements.append(sec_exp)

                formatted_list = []
                for sec_exp in tuple_elements:
                    tuple_el = nfgen.TupleElementForOutput(
                        qualifier=sub_qual, expression=sec_exp
                    )
                    formatted_list.append(tuple_el.get_string())

                expression = ", ".join(formatted_list)
                qual = nfgen.OutputProcessQualifier.tuple

            out = nfgen.ProcessOutput(
                qualifier=qual,
                name=o.id(),
                expression=expression,
                # is_optional=output_type.optional, # disable this because nextflow doesn't allow workflow to point to optional output
            )

            nf_outputs.append(out)

        return nf_outputs

    @classmethod
    def gen_output_expression(cls, o: Union[TOutput, ToolOutput]):
        """
        Based on the Janis output type, we generate string expression to represent outputs in Nextflow Process.

        :param tool:
        :type tool:
        :param o:
        :type o:
        :return:
        :rtype:
        """
        if isinstance(o, TOutput):
            output_type = o.outtype
        elif isinstance(o, ToolOutput):
            output_type = o.output_type
        else:
            raise Exception("Unknown output object")

        if isinstance(output_type, File):
            expression = f"'*{output_type.extension}'"
            qual = nfgen.OutputProcessQualifier.path
        elif isinstance(output_type, Array) and isinstance(output_type.subtype(), File):
            expression = f"'*{output_type.subtype().extension}'"
            qual = nfgen.OutputProcessQualifier.path
        else:
            qual = nfgen.OutputProcessQualifier.val
            expression = f"file(\"$workDir/{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{o.tag}\").text.replace('[', '').replace(']', '').split(', ')"

        return qual, expression

    @classmethod
    def gen_process_from_codetool(
        cls,
        tool: CommandTool,
        values: dict[str, Any],   # values fed to tool inputs (step translation)
        scope: list[str],
    ) -> nfgen.Process:
        """
        Generate a Nextflow Process object for Janis python code tool

        :param tool:
        :type tool:
        :param name:
        :type name:
        :param provided_inputs:
        :type provided_inputs:
        :return:
        :rtype:
        """
        raise NotImplementedError
        step_inputs = cls.gen_step_inval_dict(step)
        inputs: List[TInput] = []
        outputs: List[TOutput] = tool.outputs()

        # construct script
        for i in step.tool.inputs():
            if step_inputs is not None:
                optional = i.intype.optional is None or i.intype.optional is True
                if i.id() in step_inputs or i.default is not None or not optional:
                    inputs.append(i)
                # note: Filename is set to optional=True, but they have a default value that is not set to i.default
                if isinstance(i.intype, Filename):
                    inputs.append(i)
            else:
                inputs.append(i)

        script = cls.prepare_script_for_python_code_tool(step.tool, inputs)

        process_name = name or step.id()
        process = nfgen.Process(
            name=process_name,
            script=script,
            script_type=nfgen.ProcessScriptType.script,
        )

        python_file_input = nfgen.ProcessInput(
            qualifier=nfgen.InputProcessQualifier.path,
            name=settings.PYTHON_CODE_FILE_PATH_PARAM.strip("%"),
            as_param=settings.PYTHON_CODE_FILE_PATH_PARAM,
        )
        process.inputs.append(python_file_input)
        for i in inputs:
            qual = cls.get_input_qualifier_for_inptype(i.intype)
            inp = nfgen.ProcessInput(qualifier=qual, name=i.id())

            if isinstance(i.intype, File) or (
                isinstance(i.intype, Array) and isinstance(i.intype.subtype(), File)
            ):
                inp.as_param = settings.LIST_OF_FILES_PARAM

            process.inputs.append(inp)

        for o in outputs:
            qual, expression = cls.gen_output_expression(o)
            out = nfgen.ProcessOutput(
                qualifier=qual, name=o.id(), expression=expression
            )
            process.outputs.append(out)

        #resource_param_names = cls.get_resource_param_names(step.tool, name)
        #resources = cls.build_resources_input(tool, hints=None)
        resources = {}
        process.directives = cls.gen_directives_for_process(
            resources, scope
        )
        return process

    @classmethod
    def handle_process_args(
        cls,
        tool: CommandTool,
        process_inputs: Union[List[nfgen.ProcessInput], List[nfgen.WorkflowInput]],
        step_inputs: Dict[str, Any],
        workflow: Workflow,
        input_param_prefix: str = "",
        step_key_prefix: str = "",
        scatter: Optional[str] = None,
    ) -> list[str]:
        """
        generate list of strings representing each argument to a Nextflow process

        :param tool:
        :type tool:
        :param process_inputs:
        :type process_inputs:
        :param inputsdict:
        :type inputsdict:
        :param provided_inputs:
        :type provided_inputs:
        :param workflow:
        :type workflow:
        :param scatter:
        :type scatter:
        :return:
        :rtype:
        """

        # TODO scatter shouldn't be calculated here anymore. 

        args_list: list[Any] = []
        tool_inputs = {x.id(): x for x in tool.inputs()}

        # process_inputs = nextflow process {} inputs.
        # work from this perspective, so we know order matches up.
        # also know all process inputs are accounted for in step inputs. 
        for process_input in process_inputs:
            value = step_inputs[process_input.name]
            
            # elif i.name in inputsdict:
            #     p = inputsdict[i.name].default
            #     if isinstance(p, str) and "inputs." in p:
            #         p = p.replace("inputs.", f"params.{input_param_prefix}").strip("'")
            # else:
            #     p = i.as_param
            
            # Extra processing
            if settings.PYTHON_CODE_FILE_PATH_PARAM in value:
                path_to_python_code_file = posixpath.join(
                    #"$baseDir", cls.DIR_TOOLS, f"{tool.versioned_id()}.py"
                    "$baseDir", cls.DIR_TOOLS, f"{tool.id()}.py"
                )
                value = value.replace(
                    settings.PYTHON_CODE_FILE_PATH_PARAM, f'"{path_to_python_code_file}"'
                )
            # elif i.name in inputsdict:
            #     toolinput = inputsdict.get(i.name)

            #     if isinstance(toolinput, ToolInput):
            #         input_type = toolinput.input_type
            #     elif isinstance(toolinput, TInput):
            #         input_type = toolinput.intype

                # if isinstance(input_type, File):
                #     if step_input.startswith("params."):
                #         step_input = f"Channel.fromPath({step_input}).collect()"
                # elif isinstance(input_type, Array) and isinstance(
                #     input_type.subtype(), Array
                # ):
                #     if step_input.startswith("params."):
                #         step_input = f"Channel.from({step_input}).map{{ pair -> pair }}"
                # # we need to concat multiple list of channels
                # elif isinstance(input_type, Array) and isinstance(
                #     input_type.subtype(), File
                # ):
                #     if input_type.subtype().has_secondary_files():
                #         step_input = step_input.strip("][").replace(",", " + ")

            if scatter is not None and process_input.name in scatter.fields:
                tool_input = tool_inputs[process_input.name]
                step_keys = [f"{step_key_prefix}{s}" for s in list(workflow.step_nodes.keys())]
                value = cls.handle_scatter_argument(value, tool_input, step_keys)

            args_list.append(value)

        return args_list

    @classmethod
    def handle_scatter_argument(
        cls, p: str, input: ToolInput, step_keys: List[str] = []
    ):
        """

        :param p: string to represent the argument
        :type p: str
        :param step_keys: List of worfklow step names
        :type step_keys: List[str]

        :return:
        :rtype:
        """
        if isinstance(input, ToolInput):
            input_type = input.input_type
        elif isinstance(input, TInput):
            input_type = input.intype
        else:
            raise Exception("Unknown input object")

        # If it comes from our input files, it is already properly formatted
        # e.g array of array are in the correct formats
        matches = []
        pattern = r"\b(params(\..+?)+)\b"
        found = re.findall(pattern, p)
        # found is in this format
        # [('params.intervals', '.intervals'), ('step_id.out.test', '.test')]
        if found is not None:
            matches += [t[0] for t in found]

        for m in matches:
            p = p.replace(m, f"Channel.from({m}).map{{ item -> item }}")

        # Handling outputs from internal workflow steps
        matches = []
        for step_id in step_keys:
            pattern = rf"\b({step_id}(\..+?)+)\b"
            found = re.findall(pattern, p)
            # found is in this format
            # [('params.intervals', '.intervals'), ('step_id.out.test', '.test')]
            if found is not None:
                matches += [t[0] for t in found]

        for m in matches or []:
            if hasattr(input_type, 'has_secondary_files'):
                if input_type.has_secondary_files() or input_type.is_paired():
                    p = p.replace(m, f"{m}.map{{ pair -> pair }}")
            else:
                p = p.replace(m, f"{m}.flatten()")

        return p

    @classmethod
    def gen_process_workflow(
        cls, tool: Tool, process: nfgen.Process
    ) -> nfgen.Workflow:
        """
        In the main translation file, we call a Nextflow Workflow even if it is only for a tool.
        This function generates this Nextflow Workflow object.

        :param tool:
        :type tool:
        :param process:
        :type process:
        :return:
        :rtype:
        """
        body: list[str] = []

        # gather input args for the tool process call
        args_list = []
        for i in process.inputs:
            p = f"ch_{i.name}"

            # Extra processing when we need to set up the process input parameters
            if i.as_param:
                if settings.LIST_OF_FILES_PARAM in i.as_param:
                    p = i.as_param.replace(
                        settings.LIST_OF_FILES_PARAM, f"Channel.fromPath({p}).collect()"
                    )
                elif settings.LIST_OF_FILE_PAIRS_PARAM in i.as_param:
                    p = i.as_param.replace(
                        settings.LIST_OF_FILE_PAIRS_PARAM,
                        f"Channel.from({p}).map{{ pair -> pair }}",
                    )
                elif settings.PYTHON_CODE_FILE_PATH_PARAM in i.as_param:
                    path_to_python_code_file = posixpath.join(
                        #"$baseDir", cls.DIR_TOOLS, f"{tool.versioned_id()}.py"
                        "$baseDir", cls.DIR_TOOLS, f"{tool.id()}.py"
                    )
                    p = i.as_param.replace(
                        settings.PYTHON_CODE_FILE_PATH_PARAM, f'"{path_to_python_code_file}"'
                    )
            args_list.append(p)

        # gather input args for the output collection process call
        #output_args_list = [f"{process.name}.out.{o.name}" for o in process.outputs]
        
        body.append(format_process_call(process.name, args_list))
        #body.append(format_process_call(settings.FINAL_STEP_NAME, output_args_list))
        return nfgen.Workflow(name="", main=body)

    @classmethod
    def gen_step_inval_dict(
        cls,
        tool: CommandTool,
        values: dict[str, Any],
        scatter: Optional[ScatterDescription],
        inputs_replacement: str = "ch_",
        tool_id_prefix: str = "$",
    ) -> dict[str, Any]:
        """
        Generate a dictionary to represent the input values provided to a tool.
        key is the input name.
        value is the input expression.
        
        Uses StepNode.sources. Only gets step input values which are channels. 

        :param tool:
        :type tool:
        :param inputs_replacement:
        :type inputs_replacement:
        :param tool_id_prefix:
        :type tool_id_prefix:
        :return:
        :rtype:
        """

        inputs = {}
        process_ids = nfgen.utils.get_process_input_ids(tool, values)
        process_inputs = nfgen.utils.items_with_id(tool.inputs(), process_ids)
        process_inputs_names = [x.id() for x in process_inputs]



        for name in process_inputs_names:
            if name in values:
                src = values[name]
                inputs[name] = nfgen.unwrap_source(src, scatter)

                #elif isinstance(node,)

                # elif hasattr(src, "nextflow"):
                #     val = src.nextflow(
                #         var_indicator=inputs_replacement, 
                #         step_indicator=tool_id_prefix
                #     )
                # else:
                #     val = cls.unwrap_expression(
                #         src,
                #         tool=tool,
                #         inputs_dict=tool.inputs_map(), ## should this be 'values' instead?
                #         quote_string=False,
                #         var_indicator=inputs_replacement,
                #         step_indicator=tool_id_prefix,
                #     )

                # inputs[name] = val
        return inputs

    @classmethod
    def gen_wf_tool_outputs(
        cls, wf: WorkflowBase, tool_var_prefix: str = ""
    ) -> Dict[str, str]:
        """
        Generate a dictionary containing values of tool output expressions
        key is the output tag name
        value is the output expression

        :param wf:
        :type wf:
        :param tool_var_prefix:
        :type tool_var_prefix:
        :return:
        :rtype:
        """

        outputs = {}
        for o in wf.output_nodes:
            if hasattr(wf.output_nodes[o].source, "nextflow"):
                val = wf.output_nodes[o].source.nextflow(step_indicator=tool_var_prefix)
            else:
                val = str(val)

            outputs[o] = val

        return outputs

    @classmethod
    def gen_generic_functions(cls) -> List[nfgen.Function]:
        """
        Generate Groovy helper functions to run our Nextflow processes

        :return:
        :rtype:
        """

        functions = [
            nfgen.Function(
                name="apply_prefix",
                parameters=["var", "prefix", "prefix_applies_to_all_elements"],
                definition=f"""
if (prefix_applies_to_all_elements == 'True') {{
    def l = var.split(' ')
    
    prefix = prefix.toString() 
    
    return prefix + l.join(' ' + prefix)
}}
else {{
    return prefix.toString() + var.toString()
}}
""",
            ),
            nfgen.Function(
                name="optional",
                parameters=["var", "prefix", "prefix_applies_to_all_elements"],
                definition=f"""
var = var.toString()
if (var && ( var != 'None' ) && (! var.contains('{settings.NO_FILE_PATH_PREFIX}')))
{{
    return apply_prefix(var, prefix, prefix_applies_to_all_elements)
}}
else
{{
    return ''
}}
""",
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
""",
            ),
            nfgen.Function(
                name="get_primary_files",
                parameters=["var"],
                definition=f"""
def primary = []

var.eachWithIndex {{item, index -> 
    if (index % 2 == 0) {{
        primary.add(item)
    }}
}}

return primary

""",
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
        var_indicator=None,
        step_indicator=None,
        **debugkwargs,
    ): 
        return nfgen.unwrap_expression(
            value,
            input_in_selectors=cls.INPUT_IN_SELECTORS,
            quote_string=quote_string,
            tool=tool,
            for_output=for_output,
            inputs_dict=inputs_dict,
            skip_inputs_lookup=skip_inputs_lookup,
            in_shell_script=in_shell_script,
            var_indicator=var_indicator,
            step_indicator=step_indicator,
            debugkwargs=debugkwargs
        )       

    @classmethod
    def build_inputs_file(
        cls,
        tool: Workflow | CommandTool,
        recursive: bool=False,
        merge_resources: bool=False,
        hints: Optional[dict[str, Any]]=None,
        additional_inputs: Optional[dict[str, Any]]=None,
        max_cores: Optional[int]=None,
        max_mem: Optional[int]=None,
        max_duration: Optional[int]=None,
    ) -> dict[str, Any]:
        """
        Generate a dictionary containing input name and its values from users or
        its default values if inputs not provided.
        Builds using params - all params must be registered before they will appear! 
        nfgen.params.register(Tool | Workflow)

        :param tool:
        :type tool:
        :param recursive:
        :type recursive:
        :param merge_resources:
        :type merge_resources:
        :param hints:
        :type hints:
        :param additional_inputs:
        :type additional_inputs:
        :param max_cores:
        :type max_cores:
        :param max_mem:
        :type max_mem:
        :param max_duration:
        :type max_duration:
        :return:
        :rtype:
        """
        scope = []
        if additional_inputs:
            nfgen.params.register_params_for_additional_inputs(additional_inputs, scope)
        if merge_resources:
            resources_input = cls.build_resources_input(
                tool,
                hints,
                max_cores=max_cores,
                max_mem=max_mem,
                max_duration=max_duration,
            )
            nfgen.params.register_params_for_resources(resources_input)
        return nfgen.params.serialize()

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf

    @staticmethod
    def stringify_translated_tool(tool):
        return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        """
        convert dictionary inputs to string
        
        NOTE - this method does not do what it says. 

        build_inputs_file() and stringify_translated_inputs() must be 
        implemented by any subclass of TranslationBase, but these don't
        properly capture what we want to do for nextflow. 

        The fed inputs are ignored because we do not want a dict. 
        Rather, we want the registered nfgen.params. We can use their 
        properties to create a nicer nextflow.config file structure.
        
        :param inputs:
        :type inputs:
        :return:
        :rtype:
        """
        return nfgen.config.generate_config()


    @staticmethod
    def workflow_filename(workflow):
        """
        Generate the main workflow filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        #return workflow.versioned_id() + ".nf"
        return workflow.id() + ".nf"

    @staticmethod
    def inputs_filename(workflow):
        """
        Generate the input filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        #return workflow.versioned_id() + ".input.json"
        return 'nextflow.config'

    @staticmethod
    def tool_filename(tool):
        prefix = tool
        if isinstance(tool, Tool):
            #prefix = tool.versioned_id()
            prefix = tool.id()

        return prefix + ".nf"

    @staticmethod
    def resources_filename(workflow):
        """
        Generate resoureces filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        return workflow.id() + "-resources.json"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass

    @classmethod
    def prepare_script_for_python_code_tool(cls, tool: CodeTool, inputs) -> str:
        """
        Generate the content of the script section in a Nextflow process for Janis python code tool

        :param tool:
        :type tool:
        :param inputs:
        :type inputs:
        :return:
        :rtype:
        """
        PYTHON_SHEBANG = "#!/usr/bin/env python"
        #python_script_filename = f"{tool.versioned_id()}"
        python_script_filename = f"{tool.id()}"

        # TODO: handle args of type list of string (need to quote them)
        all_args = []
        for i in inputs:
            arg_value = f"${i.tag}"
            if isinstance(i.intype, Array):
                arg_value = f'"{arg_value}".split(" ")'
            elif isinstance(i.intype, File) and i.intype.has_secondary_files():
                arg_value = f'"{arg_value}".split(" ")[0]'
            elif not isinstance(i.intype, (Array, Int, Float, Double)):
                arg_value = f'"{arg_value}"'

            arg_value = f"{i.tag}={arg_value}"

            all_args.append(arg_value)

        args = ", ".join(a for a in all_args)

        script = f"""{PYTHON_SHEBANG}
from {python_script_filename} import code_block
import os
import json

result = code_block({args})

work_dir = os.getenv("PYENV_DIR")
for key in result:
    with open(os.path.join("$workDir", f"{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{{key}}"), "w") as f:
        f.write(json.dumps(result[key]))
"""

        return script


    # @classmethod
    # def prepare_tool_output(cls, tool: Tool) -> Dict[str, str]:
    #     """
    #     Generate a dictionary that contains Nextflow expressions to represent Janis outputs

    #     :param tool:
    #     :type tool:
    #     :return:
    #     :rtype:
    #     """
    #     outputs = {}
    #     for out in tool.outputs():
    #         if isinstance(out, TOutput):
    #             output_type = out.outtype
    #         elif isinstance(out, ToolOutput):
    #             output_type = out.output_type

    #         val = f"{tool.id()}{out.tag}"

    #         if (
    #             output_type.is_array()
    #             and isinstance(output_type.subtype(), File)
    #             and output_type.subtype().has_secondary_files()
    #         ):
    #             val = f"{val}.map{{ item -> item[0] }}"
    #         elif isinstance(output_type, File) and output_type.has_secondary_files():
    #             val = f"{val}[0]"

    #         outputs[out.tag] = f"${{{val}}}"

    #     return outputs

    @classmethod
    def get_input_qualifier_for_inptype(
        cls, inp_type: DataType
    ) -> nfgen.InputProcessQualifier:
        """
        Get Nextflow input qualifier based on Janis data type

        :param inp_type:
        :type inp_type:
        :return:
        :rtype:
        """

        if isinstance(inp_type, Array):
            inp_type = inp_type.fundamental_type()

        if isinstance(inp_type, (File, Directory)):
            return nfgen.InputProcessQualifier.path

        # Handle UnionType
        if inp_type.is_base_type(File) or inp_type.is_base_type(Directory):
            return nfgen.InputProcessQualifier.path

        return nfgen.InputProcessQualifier.val

    @classmethod
    def get_output_qualifier_for_outtype(
        cls,
        out_type: DataType,
    ) -> nfgen.OutputProcessQualifier:
        """
        Generate Nextflow output qualifier based on Janis output data type

        :param out_type:
        :type out_type:
        :return:
        :rtype:
        """
        if isinstance(out_type, Array):
            return nfgen.OutputProcessQualifier.tuple

        # elif isinstance(out_type, Stdout):
        #     return nfgen.OutputProcessQualifier.stdout

        elif isinstance(out_type, (File, Directory, Stdout)):
            return nfgen.OutputProcessQualifier.path

        return nfgen.OutputProcessQualifier.val

