from copy import deepcopy
import os
import re
import posixpath
from collections import defaultdict
from typing import Protocol

from typing import Tuple, Dict, List, Optional, Any
from janis_core.translations.nfgen import format_process_call

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
from janis_core.workflow.workflow import Workflow, WorkflowBase
from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.translations import nfgen
from janis_core.translations.nfgen import settings
from janis_core.translations.nfgen import Scope

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


class IGetStringMethod(Protocol):
    def get_string(self) -> str:
        ... 

def dot_to_scope_notation(dot: str) -> list[str]:
    return dot.split('.')


class NFFileRegister:
    # TODO CHECK
    """
    Stores generated nf files. 
    Organised by {scope: file}.
    Each entity (ie workflow, subworkflow, tool) only will have 1 file.
    """
    def __init__(self):
        self.files: dict[str, nfgen.NFFile] = {}

    def add(self, scope: Scope, nf_file: nfgen.NFFile) -> None:
        label = scope.to_string()
        self.files[label] = nf_file
    
    def get(self, scope: Scope) -> nfgen.NFFile:
        label = scope.to_string()
        return self.files[label]

    def get_children(self, scope: Scope, direct_only: bool=True) -> list[nfgen.NFFile]:
        # items with current scope
        child_files: list[nfgen.NFFile] = []
        
        # items with child scope
        depth = len(scope.items)
        for label, nf_file in self.files.items():
            label_split = dot_to_scope_notation(label)
            # ignore the main workflow file (it throws things off)
            if label_split != [settings.NF_MAIN_NAME]:
                # does the scope match the start of the label?
                if label_split[:depth] == scope.labels:
                    # if we only want direct children
                    if not direct_only:
                        child_files.append(nf_file)
                    elif direct_only and len(label_split) == depth + 1:
                        child_files.append(nf_file)
        return child_files
    

class NFItemRegister:
    """
    Stores generated nf items. 
    Organised by {scope: [items]}.
    Some entities (Tools / Workflows) may have multiple items. 
    """
    def __init__(self):
        self.items: dict[str, list[IGetStringMethod]] = defaultdict(list)

    def add(self, scope: Scope, nf_item: IGetStringMethod) -> None:
        label = scope.to_string()
        self.items[label].append(nf_item)

    def get(self, scope: Scope) -> list[IGetStringMethod]:
        # items with current scope
        label = scope.to_string()
        nf_items = deepcopy(self.items[label]) # ???
        return nf_items
    

class NextflowTranslator(TranslatorBase):
    DIR_TOOLS: str = '' # DO NOT ALTER
    SUBDIRS_TO_CREATE: list[str] = [
        settings.PROCESS_OUTDIR,
        settings.SUBWORKFLOW_OUTDIR,
        settings.CODE_FILES_OUTDIR,
    ]

    file_register: NFFileRegister = NFFileRegister()
    item_register: NFItemRegister = NFItemRegister()

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

        # set class variables to avoid passing junk params
        settings.BASE_OUTDIR = cls.basedir
        settings.WITH_CONTAINER = with_container
        settings.ALLOW_EMPTY_CONTAINER = allow_empty_container
        settings.CONTAINER_OVERRIDE = container_override
        settings.WITH_RESOURCE_OVERRIDES = with_resource_overrides
        settings.RENDER_COMMENTS = render_comments

        # blank scope - main wf has not parent
        scope = nfgen.Scope()

        # register params and channels for workflow inputs
        nfgen.params.clear()
        nfgen.channels.clear()
        nfgen.register_params_channels(jworkflow, scope)

        # main logic
        cls.update_files('', scope, jworkflow, sources={})

        # get the main wf file and all sub files
        main_file = cls.file_register.get(scope)  # main nf workflow usually
        sub_files = cls.file_register.get_children(scope, direct_only=False)

        # return format (gen str for each file)
        main_file_str = main_file.get_string()
        sub_files_str = {sub_file.path: sub_file.get_string() for sub_file in sub_files}
        return (main_file_str, sub_files_str)


    @classmethod
    def update_files(
        cls, 
        identifier: str,
        scope: Scope,
        tool: Workflow | CommandTool | PythonTool,
        sources: dict[str, Any],
        scatter: Optional[ScatterDescription]=None,
        ) -> None:
        
        """
        1. generate all nextflow items (currently: nfgen.Process, nfgen.Workflow, nfgen.ChannelOperation).
        2. write any nextflow items that are processes or workflows to file.
        3. return the nextflow items (so a workflow in scope above can generate process / workflow calls).
        """
        subtype: str = nfgen.naming.get_construct_name(tool, scope)

        # groovy code for plumbing
        # TODO: reimplement
        # if scatter and scatter.method == ScatterMethod.cross:
        #     operation_item = nfgen.channels.gen_scatter_cross_operation(sources, scatter)
        #     cls.item_register.add(scope, operation_item)

        # command tool
        if isinstance(tool, CommandTool):
            # item
            process_item = nfgen.process.gen_process_from_cmdtool(tool, sources, scope)
            process_item = cls.handle_container(tool, process_item)
            cls.item_register.add(scope, process_item)
            # file
            process_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=[], 
                items=cls.item_register.get(scope),
            )
            cls.file_register.add(scope, process_file)

        # python tool
        elif isinstance(tool, PythonTool):
            # item
            process_item = nfgen.process.gen_process_from_codetool(tool, sources, scope)
            process_item = cls.handle_container(tool, process_item)
            cls.item_register.add(scope, process_item)
            # file
            process_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=[], 
                items=cls.item_register.get(scope),
            )
            cls.file_register.add(scope, process_file)

        # workflow
        elif isinstance(tool, WorkflowBase):
            
            # handle sub elements (tool / workflow calls)
            wf = tool
            for substep in wf.step_nodes.values():
                current_scope = deepcopy(scope)
                current_scope.update(substep)
                cls.update_files(
                    identifier=current_scope.labels[-1],
                    scope=current_scope,
                    tool=substep.tool,
                    sources=substep.sources,
                    scatter=substep.scatter
                )

            # item: channels item (if main workflow object)
            if scope.labels == [settings.NF_MAIN_NAME]:
                channels_item = nfgen.channels.channels.channel_register # bad naming
                if len(channels_item.ordered_channels) > 0:
                    cls.item_register.add(scope, channels_item)

            # item: workflow body
            workflow_item = nfgen.workflow.gen_workflow(
                name=identifier,
                scope=scope,
                sources=sources,
                wf=wf,
                item_register=cls.item_register
            )
            cls.item_register.add(scope, workflow_item)

            # imports
            imports: list[nfgen.Import] = []
            for nf_file in cls.file_register.get_children(scope):
                # get the relative import path. this is ugly last-min code. 
                # (main wf)
                if len(scope.labels) == 1: 
                    source = f'./{nf_file.path}'
                # (subwf)
                elif settings.SUBWORKFLOW_OUTDIR in nf_file.path:
                    folder = os.path.split(nf_file.path)[-1]
                    source = f'./{folder}'
                elif settings.PROCESS_OUTDIR in nf_file.path:
                    source = f'../{nf_file.path}'
                else:
                    raise NotImplementedError
                
                nf_item = nfgen.ImportItem(name=nf_file.name)
                nf_import = nfgen.Import(items=[nf_item], source=source)
                imports.append(nf_import)

            # file: workflow file            
            workflow_file = nfgen.NFFile(
                name=identifier,  # TODO here name clash checking
                subtype=subtype,
                imports=imports, 
                items=cls.item_register.get(scope),
            )
            cls.file_register.add(scope, workflow_file)

        else:
            raise Exception(
                f"Nextflow translation for this {tool.__class__} is not yet supported"
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
        raise NotImplementedError
        settings.MODE = 'tool'
        file_items: list[IGetStringMethod] = []
        scope: Scope = nfgen.Scope()
        values: dict[str, Any] = {}

        main_p = nfgen.process.gen_process_from_cmdtool(tool, values, scope)
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
        if settings.WITH_CONTAINER:
            container = (
                NextflowTranslator.get_container_override_for_tool(
                    tool, settings.CONTAINER_OVERRIDE
                )
                or tool.container()
            )

        if container is not None:
            container_expr = nfgen.unwrap_expression(container, tool=tool)
            directive = nfgen.ContainerDirective(container_expr)
            process.directives.append(directive)

        elif not settings.ALLOW_EMPTY_CONTAINER:
            raise Exception(
                f"The tool '{tool.id()}' did not have a container and no container override was specified. "
                f"Although not recommended, Janis can export empty docker containers with the parameter "
                f"'allow_empty_container=True' or --allow-empty-container"
            )

        return process

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
            # helpers["__init__.py"] = ""
            #helpers[f"{tool.versioned_id()}.py"] = cls.gen_python_script(tool)
            subdir = settings.CODE_FILES_OUTDIR
            filename = f'{tool.id()}.py'
            filepath = os.path.join(subdir, filename)
            helpers[filepath] = tool.prepared_script(SupportedTranslation.NEXTFLOW)
            return helpers

        elif isinstance(tool, WorkflowBase):
            for step_id in tool.step_nodes:
                step_tool = tool.step_nodes[step_id].tool
                helpers = cls.gen_python_code_files(step_tool, helpers)
        return helpers

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
        process_ids = nfgen.process.inputs.get_process_inputs(sources)
        scatter_method = scatter.method if scatter else None

        for name in process_ids:
            if name in sources:
                src = sources[name]
                scatter_target = True if scatter and name in scatter.fields else False

                inputs[name] = nfgen.unwrap_expression(
                    val=src, 
                    tool=tool,
                    sources=sources,
                    scatter_target=scatter_target,
                    scatter_method=scatter_method,
                )
        return inputs

    @classmethod
    def unwrap_expression(
        cls,
        value,
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
            tool=tool,
            # inputs_dict=inputs_dict,
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
        scope: Scope = nfgen.Scope()
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




    # DEPRECATED

    # @classmethod
    # def collect_workflow_imports(
    #     cls, files: dict[str, nfgen.NFFile]
    # ) -> list[nfgen.Import]:
    #     """
    #     Generate a list of Nextflow Import objects for each of the steps in a workflow.

    #     :param files:
    #     :type files:
    #     :return:
    #     :rtype:
    #     """
    #     imports: list[nfgen.Import] = []
        
    #     for filename, nf_file in files.items():
    #         file_imports = []
    #         for process in nf_file.items:
    #             i_item = nfgen.ImportItem(name=process.name)
    #             file_imports.append(i_item)

    #         imp = nfgen.Import(
    #             file_imports, os.path.join(".", settings.PROCESS_OUTDIR, filename)
    #         )
    #         imports.append(imp)

    #     return imports


    # @classmethod
    # def handle_scatter_argument(
    #     cls, p: str, input: ToolInput, step_keys: List[str] = []
    # ):
    #     """

    #     :param p: string to represent the argument
    #     :type p: str
    #     :param step_keys: List of worfklow step names
    #     :type step_keys: List[str]

    #     :return:
    #     :rtype:
    #     """
    #     if isinstance(input, ToolInput):
    #         input_type = input.input_type
    #     elif isinstance(input, TInput):
    #         input_type = input.intype
    #     else:
    #         raise Exception("Unknown input object")

    #     # If it comes from our input files, it is already properly formatted
    #     # e.g array of array are in the correct formats
    #     matches = []
    #     pattern = r"\b(params(\..+?)+)\b"
    #     found = re.findall(pattern, p)
    #     # found is in this format
    #     # [('params.intervals', '.intervals'), ('step_id.out.test', '.test')]
    #     if found is not None:
    #         matches += [t[0] for t in found]

    #     for m in matches:
    #         p = p.replace(m, f"Channel.from({m}).map{{ item -> item }}")

    #     # Handling outputs from internal workflow steps
    #     matches = []
    #     for step_id in step_keys:
    #         pattern = rf"\b({step_id}(\..+?)+)\b"
    #         found = re.findall(pattern, p)
    #         # found is in this format
    #         # [('params.intervals', '.intervals'), ('step_id.out.test', '.test')]
    #         if found is not None:
    #             matches += [t[0] for t in found]

    #     for m in matches or []:
    #         if hasattr(input_type, 'has_secondary_files'):
    #             if input_type.has_secondary_files() or input_type.is_paired():
    #                 p = p.replace(m, f"{m}.map{{ pair -> pair }}")
    #         else:
    #             p = p.replace(m, f"{m}.flatten()")

    #     return p

    # @classmethod
    # def gen_process_workflow(
    #     cls, tool: Tool, process: nfgen.Process
    # ) -> nfgen.Workflow:
    #     """
    #     In the main translation file, we call a Nextflow Workflow even if it is only for a tool.
    #     This function generates this Nextflow Workflow object.

    #     :param tool:
    #     :type tool:
    #     :param process:
    #     :type process:
    #     :return:
    #     :rtype:
    #     """
    #     main: list[str] = []
    #     name = ""

    #     # gather input args for the tool process call
    #     args_list = []
    #     for i in process.inputs:
    #         p = f"ch_{i.name}"

    #         # Extra processing when we need to set up the process input parameters
    #         if i.as_param:
    #             if settings.LIST_OF_FILES_PARAM in i.as_param:
    #                 p = i.as_param.replace(
    #                     settings.LIST_OF_FILES_PARAM, f"Channel.fromPath({p}).collect()"
    #                 )
    #             elif settings.LIST_OF_FILE_PAIRS_PARAM in i.as_param:
    #                 p = i.as_param.replace(
    #                     settings.LIST_OF_FILE_PAIRS_PARAM,
    #                     f"Channel.from({p}).map{{ pair -> pair }}",
    #                 )
    #             elif settings.PYTHON_CODE_FILE_SYMBOL in i.as_param:
    #                 path_to_python_code_file = posixpath.join(
    #                     #"$baseDir", settings.PROCESS_OUTDIR, f"{tool.versioned_id()}.py"
    #                     "$baseDir", settings.PROCESS_OUTDIR, f"{tool.id()}.py"
    #                 )
    #                 p = i.as_param.replace(
    #                     settings.PYTHON_CODE_FILE_SYMBOL, f'"{path_to_python_code_file}"'
    #                 )
    #         args_list.append(p)

    #     # gather input args for the output collection process call
    #     #output_args_list = [f"{process.name}.out.{o.name}" for o in process.outputs]
        
    #     main.append(format_process_call(process.name, args_list))
    #     #body.append(format_process_call(settings.FINAL_STEP_NAME, output_args_list))
    #     raise NotImplementedError
    #     return nfgen.Workflow(name, main)


    # @classmethod
    # def gen_wf_tool_outputs(
    #     cls, wf: WorkflowBase, tool_var_prefix: str = ""
    # ) -> Dict[str, str]:
    #     """
    #     Generate a dictionary containing values of tool output expressions
    #     key is the output tag name
    #     value is the output expression

    #     :param wf:
    #     :type wf:
    #     :param tool_var_prefix:
    #     :type tool_var_prefix:
    #     :return:
    #     :rtype:
    #     """
    #     outputs = {}
    #     for o in wf.output_nodes:
    #         if hasattr(wf.output_nodes[o].source, "nextflow"):
    #             val = wf.output_nodes[o].source.to_nextflow(step_indicator=tool_var_prefix)
    #         else:
    #             val = str(val)
    #         outputs[o] = val
    #     return outputs