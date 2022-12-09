from copy import deepcopy
import os
import re
import posixpath
from collections import defaultdict

from typing import Tuple, Dict, List, Optional, Any
from janis_core.translations.nfgen import format_process_call
from janis_core.operators.selectors import Selector
from janis_core.types import (
    Array,
    Int,
    Float,
    Double,
)

from janis_core.utils.scatter import ScatterDescription, ScatterMethod

from janis_core.tool.commandtool import (
    CommandTool,
    ToolInput,
    Tool,
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

# - GH Dec 2022


class NextflowTranslator(TranslatorBase):
    DIR_TOOLS: str = "" # dont change
    SUBDIRS_TO_CREATE: list[str] = [
        'modules',
        'subworkflows',
    ]
    nf_files: dict[str, nfgen.NFFile] = {}
    nf_items: dict[str, list[nfgen.NFBase]] = defaultdict(list)

    def __init__(self):
        super().__init__(name="nextflow")
        # self.nf_items: dict[str, list[nfgen.NFBase]] = defaultdict(list)
        # self.files: dict[str, str] = {}

    @classmethod
    def scope_to_dot_notation(cls, scope: list[str]) -> str:
        if not scope:
            label = ''
        else:
            label = '.'.join(scope)
        return label

    @classmethod
    def add_to_nf_files(cls, scope: list[str], nf_file: nfgen.NFFile) -> None:
        # TODO CHECK
        label = cls.scope_to_dot_notation(scope)
        cls.nf_files[label] = nf_file
    
    @classmethod
    def add_to_nf_items(cls, scope: list[str], nf_item: nfgen.NFBase) -> None:
        # TODO CHECK
        label = cls.scope_to_dot_notation(scope)
        cls.nf_items[label].append(nf_item)
    
    @classmethod
    def get_nf_file(cls, scope: list[str]) -> nfgen.NFFile:
        label = cls.scope_to_dot_notation(scope)
        return cls.nf_files[label]
    
    @classmethod
    def get_nf_items(cls, scope: list[str], include_direct_children: bool=False) -> list[nfgen.NFBase]:
        label = cls.scope_to_dot_notation(scope)
        
        # items with current scope
        nf_items = deepcopy(cls.nf_items[label])
        
        # items with child scope
        if include_direct_children:
            levels = len(scope)
            for label, items in cls.nf_items.items():
                label_split = label.split('.')
                # ignore the main workflow file (it throws things off)
                if label_split != ['']:
                    # does the scope match the start of the label?
                    if label_split[:levels] == scope:
                        # is the label only 1 level below? 
                        if len(label_split) == levels + 1:
                            nf_items += items
        return nf_items

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

        # avoiding passing junk params
        cls.with_container = with_container
        cls.with_resource_overrides = with_resource_overrides
        cls.allow_empty_container = allow_empty_container
        cls.container_override = container_override
        cls.render_comments = render_comments

        scope: list[str] = []

        # register params and channels for workflow inputs
        nfgen.register_params_channels(jworkflow, scope)

        cls.update_files(scope, jworkflow, sources={})

        # # parse each step to a NFFile
        # for step in jworkflow.step_nodes.values():
        #     current_scope = deepcopy(scope)
        #     current_scope.append(step.id())       
        #     nf_items = cls.gen_items_for_step(step, scope=current_scope)
        #     nf_file = nfgen.NFFile(
        #         subtype=nfgen.get_construct_name(step.tool),
        #         imports=[], 
        #         items=[x[1] for x in nf_items],
        #         name=step.id()  # TODO here name clash checking
        #     )
        #     self.add_to_nf_files(scope, nf_file)

        #     files[nf_file.path] = nf_file

        # # handle imports
        # imports = cls.collect_workflow_imports(files)
        
        # # create object & NFFile for workflow
        # workflow = cls.gen_workflow(
        #     wf=jworkflow, 
        #     nf_items=nf_items,
        # )
        # channels = nfgen.channels.channel_register
        # nf_file = nfgen.NFFile(subtype='', imports=imports, items=[channels, workflow])

        # generate strings for each file
        tool_scripts: dict[str, str] = {path: nffile.get_string() for path, nffile in files.items()}
        return (nf_file.get_string(), tool_scripts)


    @classmethod
    def update_files(
        cls, 
        scope: list[str],
        tool: Workflow | CommandTool | PythonTool,
        sources: dict[str, Any],
        scatter: Optional[ScatterDescription]=None,
        ) -> list[nfgen.NFBase]:
        
        """
        1. generate all nextflow items (currently: nfgen.Process, nfgen.Workflow, nfgen.ChannelOperation).
        2. write any nextflow items that are processes or workflows to file.
        3. return the nextflow items (so a workflow in scope above can generate process / workflow calls).
        """
        identifier: str = scope[-1] if scope else '__main_file__'
        subtype: str = nfgen.get_construct_name(tool, scope)

        # any groovy code
        if scatter and scatter.method == ScatterMethod.cross:
            operation_item = nfgen.channels.gen_scatter_cross_operation(sources, scatter)
            cls.add_to_nf_items(scope, operation_item)

        # command tool
        if isinstance(tool, CommandTool):
            # item
            process_item = cls.gen_process_from_cmdtool(tool, sources, scope)
            process_item = cls.handle_container(tool, process_item)
            cls.add_to_nf_items(scope, process_item)
            # file
            process_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=[], 
                items=cls.get_nf_items(scope),
            )
            cls.add_to_nf_files(scope, process_file)
            print()

        # python tool
        elif isinstance(tool, PythonTool):
            # item
            process_item = cls.gen_process_from_codetool(tool, sources, scope)
            process_item = cls.handle_container(tool, process_item)
            cls.add_to_nf_items(scope, process_item)
            # file
            process_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=[], 
                items=cls.get_nf_items(scope),
            )
            cls.add_to_nf_files(scope, process_file)
            print()

        # workflow
        elif isinstance(tool, WorkflowBase):
            
            # handle sub elements (tool / workflow calls)
            wf = tool
            for substep in wf.step_nodes.values():
                current_scope = deepcopy(scope)
                current_scope.append(substep.id())
                cls.update_files(
                    scope=current_scope,
                    tool=substep.tool,
                    sources=substep.sources,
                    scatter=substep.scatter
                )

            # item: channels item (if main workflow object)
            if identifier == '__main_file__':
                channels_item = nfgen.channels.channel_register
                cls.add_to_nf_items(scope, channels_item)

            # item: workflow body
            workflow_item = cls.gen_workflow(
                name=identifier,
                scope=scope,
                wf=wf
            )
            cls.add_to_nf_items(scope, workflow_item)

            # imports
            imports: list[nfgen.Import] = []
            for sub_file in sub_nf_files:
                nf_item = nfgen.ImportItem(name=sub_file.name)
                nf_import = nfgen.Import(items=[nf_item], source=sub_file.path)
                imports.append(nf_import)

            # file: workflow file            
            workflow_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=imports, 
                items=sub_nf_items,
            )
            cls.add_to_nf_files(scope, workflow_file)
            print()

        else:
            raise Exception(
                f"Nextflow translation for this {tool.__class__} is not yet supported"
            )


    # @classmethod
    # def gen_items_for_step(
    #     cls,
    #     step: StepNode,
    #     scope: list[str],
    # ) -> list[Tuple[str, nfgen.NFBase]]:
    #     """
    #     For each of the workflow step, we need to generate a Nextflow subworkflow or process object

    #     :param tool:
    #     :type tool:
    #     :param step_id:
    #     :type step_id:
    #     :param step_keys:
    #     :type step_keys:
    #     :param with_container:
    #     :type with_container:
    #     :param with_resource_overrides:
    #     :type with_resource_overrides:
    #     :param allow_empty_container:
    #     :type allow_empty_container:
    #     :param container_override:
    #     :type container_override:
    #     :return:
    #     :rtype:
    #     """

    #     nf_items: dict[str, list[nfgen.NFBase]] = defaultdict(list)
    #     if step.scatter and step.scatter.method == ScatterMethod.cross:
    #         operation = nfgen.channels.gen_scatter_cross_operation(step.sources, step.scatter)

    #         nf_items.append((step.id(), operation))

    #     if isinstance(step.tool, CommandTool):
    #         process = cls.gen_process_from_cmdtool(step.tool, step.sources, scope)
    #         process = cls.handle_container(step.tool, process)
    #         nf_items.append((step.id(), process))
    #         return nf_items

    #     elif isinstance(step.tool, PythonTool):
    #         process = cls.gen_process_from_codetool(step.tool, step.sources, scope)
    #         process = cls.handle_container(step.tool, process)
    #         nf_items.append((step.id(), process))
    #         return nf_items

    #     elif isinstance(step.tool, WorkflowBase):
    #         subworkflow = step.tool
    #         sub_nf_items: list[Tuple[str, nfgen.NFBase]] = []
    #         for substep in subworkflow.step_nodes.values():
    #             current_scope = deepcopy(scope)
    #             current_scope.append(substep.id())
    #             sub_nf_items += cls.gen_items_for_step(substep, current_scope)
    #         workflow_item = cls.gen_workflow(
    #             # TODO IMPORTS HERE PLS
    #             # scope=scope,
    #             wf=subworkflow, 
    #             nf_items=sub_nf_items, 
    #             name=step.id(), 
    #             is_subworkflow=True
    #         )
    #         nf_items.append((step.id(), workflow_item))
    #         nf_items += sub_nf_items
    #         return nf_items

    #     elif isinstance(step.tool, CodeTool):
    #         raise Exception("Only PythonTool code tool is supported for the moment")
        
    #     else:
    #         raise Exception(
    #             f"Nextflow translation for this {step.tool.__class__} is not yet supported"
    #         )

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
        main_p = cls.handle_container(tool, main_p)
        # output_p_inp, output_p_out = cls.prepare_output_process_params_for_tool(tool, main_p)
        # output_p = cls.gen_output_process(inputs=output_p_inp, tool_outputs=output_p_out)
        workflow = cls.gen_process_workflow(tool, process=main_p)

        file_items.append(param_block)
        #file_items.append(channel_block)
        file_items.append(main_p)
        #file_items.append(output_p)
        file_items.append(workflow)

        nf_file = nfgen.NFFile(
            subtype='',
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
            process = cls.gen_process_from_codetool(tool=tool, sources={}, scope=[])
            process = cls.handle_container(tool, process)

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
            nf_file = nfgen.NFFile(subtype='', imports=imports, items=items)

            return nf_file.get_string()
        else:
            raise Exception("Only PythonTool code tool is supported for the moment.")
        
    @classmethod
    def gen_workflow(
        cls,
        name: str,
        scope: list[str],
        wf: WorkflowBase,
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
        # TODO IMPORTS HERE PLS
        is_subworkflow = True if scope else False

        take: list[nfgen.WorkflowTake] = []
        emit: list[nfgen.WorkflowEmit] = []
        main: list[str] = []

        if is_subworkflow:
            # take
            for inp in wf.input_nodes.values():
                if nfgen.channels.exists(inp.uuid):
                    ch = nfgen.channels.get(inp.uuid)
                    take.append(nfgen.WorkflowTake(ch.name))
            
            # emit
            emit: list[nfgen.WorkflowEmit] = []
            for out in wf.output_nodes.values():
                outname = out.id()
                expression = nfgen.unwrap_expression(val=out.source)
                emit.append(nfgen.WorkflowEmit(outname, expression))

        # main
        nf_items = cls.get_nf_items(scope, include_direct_children=True)

        for step_id, nf_item in nf_items:
            step = wf.step_nodes[step_id]

            if isinstance(nf_item, nfgen.ChannelOperation):
                main.append(nf_item.get_string())
                continue
            
            elif isinstance(nf_item, nfgen.Process):
                entity_name = nf_item.name
                entity_inputs = nf_item.inputs

            elif isinstance(nf_item, nfgen.Workflow):
                entity_name = nf_item.name
                entity_inputs = nf_item.take

            else:
                raise NotImplementedError
            
            args = nfgen.get_args(step.tool, step.sources, step.scatter)
            main.append(format_process_call(entity_name, args))

        return nfgen.Workflow(name, main, take, emit, is_subworkflow)
        
    @classmethod
    def handle_container(
        cls,
        tool: Tool,
        process: nfgen.Process,
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
        if cls.with_container:
            container = (
                NextflowTranslator.get_container_override_for_tool(
                    tool, cls.container_override
                )
                or tool.container()
            )

        if container is not None:
            container_expr = nfgen.unwrap_expression(container, tool=tool, quote_string=False)
            directive = nfgen.ContainerDirective(container_expr)
            process.directives.append(directive)

        elif not cls.allow_empty_container:
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
        sources: dict[str, Any],   # values fed to tool inputs (step translation)
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
        # name
        process_name = scope[-1] if scope else tool.id()

        # directives
        resources = {}
        process_directives = cls.gen_directives_for_process(
            tool, resources, scope
        )

        # inputs
        process_inputs: list[nfgen.ProcessInput] = []
        tinput_ids = nfgen.process.get_process_inputs(sources)
        tinputs = nfgen.nfgen_utils.items_with_id(tool.inputs(), tinput_ids)
        tinputs = nfgen.ordering.order_janis_process_inputs(tinputs)
        for i in tinputs:
            process_inputs += nfgen.process.create_inputs(i)

        # outputs
        process_outputs: list[nfgen.ProcessOutput] = []
        for out in tool.outputs():
            process_outputs += nfgen.process.create_outputs(out, tool, sources)

        # script
        pre_script, script = nfgen.process.gen_script_for_cmdtool(
            tool=tool,
            scope=scope,
            sources=sources,
            stdout_filename=settings.TOOL_STDOUT_FILENAME
        )
        
        # process
        process = nfgen.Process(
            name=process_name,
            pre_script=pre_script,
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            inputs=process_inputs,
            outputs=process_outputs,
            directives=process_directives
        )

        return process

    @classmethod
    def gen_directives_for_process(
        cls, tool: CommandTool, resources: dict[str, Any], scope: Optional[list[str]]=None
    ) -> list[nfgen.ProcessDirective]:
        scope = scope if scope else []

        nf_directives: dict[str, nfgen.ProcessDirective] = {}
        nf_directives['publishDir'] = nfgen.PublishDirDirective(scope)
        nf_directives['debug'] = nfgen.DebugDirective(debug='true')

        # Add directives from input resources

        for res, val in resources.items():
            if res.endswith("runtime_cpu"):
                nf_directives['cpus'] = nfgen.CpusDirective(scope, res, val)
            elif res.endswith("runtime_memory"):
                nf_directives['memory'] = nfgen.MemoryDirective(scope, res, val)
            elif res.endswith("runtime_seconds"):
                nf_directives['time'] = nfgen.TimeDirective(scope, res, val)
            elif res.endswith("runtime_disk"):
                nf_directives['disk'] = nfgen.DiskDirective(scope, res, val)
        
        # Add directives from tool resources
        if 'cpus' not in nf_directives and tool.cpus({}) is not None:           
            nf_directives['cpus'] = nfgen.CpusDirective(scope, 'cpus', tool.cpus({}))
        if 'memory' not in nf_directives and tool.memory({}) is not None:
            nf_directives['memory'] = nfgen.MemoryDirective(scope, 'memory', tool.memory({}))
        if 'disk' not in nf_directives and tool.disk({}) is not None:
            nf_directives['disk'] = nfgen.DiskDirective(scope, 'disk', tool.disk({}))
        if 'time' not in nf_directives and tool.time({}) is not None:
            nf_directives['time'] = nfgen.TimeDirective(scope, 'time', tool.time({}))
        
        final_directives: list[nfgen.ProcessDirective] = []
        for direc in nf_directives.values():
            if hasattr(direc, 'default') and isinstance(direc.default, Selector):
                continue
            else:
                final_directives.append(direc)
        return final_directives

    @classmethod
    def gen_process_from_codetool(
        cls,
        tool: PythonTool,
        sources: dict[str, Any],   # values fed to tool inputs (step translation)
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
        # name
        process_name = scope[-1] if scope else tool.id()

        # directives
        resources = {}
        process_directives = cls.gen_directives_for_process(tool, resources, scope)
        
        # inputs
        process_inputs: list[nfgen.ProcessInput] = []
        
        # inputs: python script
        python_file_input = nfgen.process.PathProcessInput(
            name=settings.PYTHON_CODE_FILE_PATH_PARAM.strip("%"),
        )
        process_inputs.append(python_file_input)
        
        # inputs: tool inputs
        tinput_ids = nfgen.process.get_process_inputs(sources)
        tinputs = nfgen.nfgen_utils.items_with_id(tool.inputs(), tinput_ids)
        tinputs = nfgen.ordering.order_janis_process_inputs(tinputs)
        for i in tinputs:
            process_inputs += nfgen.process.create_inputs(i)

        # outputs
        process_outputs: list[nfgen.ProcessOutput] = []
        for out in tool.outputs():
            process_outputs += nfgen.process.create_outputs(out, tool, sources)

        # script
        script = cls.prepare_script_for_python_code_tool(tool, sources=sources)

        # process
        process = nfgen.Process(
            name=process_name,
            script=script,
            script_type=nfgen.ProcessScriptType.script,
            inputs=process_inputs,
            outputs=process_outputs,
            directives=process_directives
        )

        return process

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
        main: list[str] = []
        name = ""

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
        
        main.append(format_process_call(process.name, args_list))
        #body.append(format_process_call(settings.FINAL_STEP_NAME, output_args_list))
        raise NotImplementedError
        return nfgen.Workflow(name, main)

    @classmethod
    def gen_step_inval_dict(
        cls,
        tool: CommandTool,
        sources: dict[str, Any],
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
        process_ids = nfgen.process.get_process_inputs(sources)

        for name in process_ids:
            if name in sources:
                src = sources[name]
                scatter_target = True if name in scatter.fields else False

                inputs[name] = nfgen.unwrap_expression(
                    val=src, 
                    tool=tool,
                    sources=sources,
                    scatter_target=scatter_target,
                    scatter_method=scatter.method,
                )
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
                val = wf.output_nodes[o].source.to_nextflow(step_indicator=tool_var_prefix)
            else:
                val = str(val)
            outputs[o] = val
        return outputs

    @classmethod
    def unwrap_expression(
        cls,
        value,
        quote_string=True,
        tool=None,
        # inputs_dict=None,
        skip_inputs_lookup=False,
        in_shell_script=False,
        var_indicator=None,
        step_indicator=None,
        **debugkwargs,
    ): 
        return nfgen.unwrap_expression(
            val=value,
            # input_in_selectors=cls.INPUT_IN_SELECTORS,
            quote_string=quote_string,
            tool=tool,
            # inputs_dict=inputs_dict,
            skip_inputs_lookup=skip_inputs_lookup,
            in_shell_script=in_shell_script,
            # var_indicator=var_indicator,
            # step_indicator=step_indicator,
            # debugkwargs=debugkwargs
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
    def workflow_filename(workflow: Workflow) -> str:
        """
        Generate the main workflow filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        #return workflow.versioned_id() + ".nf"
        # return workflow.id() + ".nf"
        return 'main.nf'

    @staticmethod
    def inputs_filename(workflow: Workflow) -> str:
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
    def tool_filename(tool: str | Tool) -> str:
        prefix: str = '' 
        if isinstance(tool, Tool):
            #prefix = tool.versioned_id()
            prefix = tool.id()
        else:
            prefix = tool

        return prefix + ".nf"

    @staticmethod
    def resources_filename(workflow: Workflow) -> str:
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
    def prepare_script_for_python_code_tool(cls, tool: PythonTool, sources: dict[str, Any]) -> str:
        """
        Generate the content of the script section in a Nextflow process for Janis python code tool

        :param tool:
        :type tool:
        :param inputs:
        :type inputs:
        :return:
        :rtype:
        """
        python_script_filename=tool.id()

        # TODO: handle args of type list of string (need to quote them)
        args: list[str] = []
        process_inputs = nfgen.process.get_process_inputs(sources)
        param_inputs = nfgen.process.get_param_inputs(sources)
        
        for inp in tool.inputs():
            if inp.id() in process_inputs or inp.id() in param_inputs:
                src = nfgen.process.get_src(inp, process_inputs, param_inputs, sources)

                value = f"${{{src}}}"
                if isinstance(inp.intype, Array):
                    value = f'"{value}".split(" ")'
                elif not isinstance(inp.intype, (Array, Int, Float, Double)):
                    value = f'"{value}"'
                arg = f"{inp.tag}={value}"
                args.append(arg)  
            elif inp.default is not None:
                arg = f"{inp.tag}={inp.default}"
                args.append(arg)  

        args_str = ", ".join(a for a in args)
        script = f"""\
{settings.PYTHON_SHEBANG}

from {python_script_filename} import code_block
import os
import json

result = code_block({args_str})

work_dir = os.getenv("PYENV_DIR")
for key in result:
    with open(os.path.join("${{task.workDir}}", f"{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{{key}}"), "w") as f:
        f.write(json.dumps(result[key]))
"""
        return script
