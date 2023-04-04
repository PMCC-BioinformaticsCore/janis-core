from copy import deepcopy
import os
from collections import defaultdict
from typing import Tuple, Optional, Any, Protocol

from janis_core.utils.scatter import ScatterDescription
from janis_core.tool.commandtool import (
    CommandTool,
    Tool,
)

from janis_core.code.codetool import CodeTool
from janis_core.code.pythontool import PythonTool
from janis_core.workflow.workflow import Workflow, WorkflowBase
from janis_core.translations.translationbase import TranslatorBase
from janis_core.translation_deps.supportedtranslations import SupportedTranslation
from janis_core import InputSelector, File, Directory
from janis_core import settings

from .model.files.files import NFFile, NFImport, NFImportItem, NFChannelsBlock, NFFunctionsBlock, NFImportsBlock
from .model.process import NFContainerDirective
from .model.process import NFProcess
from .model.workflow import NFWorkflow

from . import channels
from . import params
from . import parsing
from . import naming
from . import preprocessing

from .scope import Scope
from .unwrap import unwrap_expression


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
    """
    Stores generated nf files. 
    Organised by {scope: file}.
    Each entity (ie workflow, subworkflow, tool) only belongs to 1 file.
    Some entities (Tools / Workflows) may have multiple items (imports, function definitions, process / workflow body etc). 
    """
    def __init__(self):
        self.files: dict[str, NFFile] = {}

    def add(self, scope: Scope, nf_file: NFFile) -> None:
        label = scope.to_string()
        self.files[label] = nf_file
    
    def get(self, scope: Scope) -> NFFile:
        label = scope.to_string()
        return self.files[label]

    def get_children(self, scope: Scope, direct_only: bool=True) -> list[NFFile]:
        # items with current scope
        child_files: list[NFFile] = []
        
        # items with child scope
        depth = len(scope.items)
        for label, nf_file in self.files.items():
            label_split = dot_to_scope_notation(label)
            # ignore the main workflow file (it throws things off)
            if label_split != [settings.translate.nextflow.NF_MAIN_NAME]:
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
    Some entities (Tools / Workflows) may have multiple items (imports, function definitions, process / workflow body etc). 
    """
    def __init__(self):
        self.items: dict[str, list[IGetStringMethod]] = defaultdict(list)

    def add(self, scope: Scope, nf_item: IGetStringMethod) -> None:
        label = scope.to_string()
        self.items[label].append(nf_item)

    def get(self, scope: Scope) -> list[IGetStringMethod]:
        # items with current scope
        label = scope.to_string()
        nf_items = deepcopy(self.items[label]) # ??? no idea why i did this but too scared to change
        return nf_items
    

class NextflowTranslator(TranslatorBase):
    DIR_TOOLS: str = '' # DO NOT ALTER
    DIR_FILES: str = '' # DO NOT ALTER
    SUBDIRS_TO_CREATE: list[str] = [
        settings.translate.nextflow.PROCESS_OUTDIR,
        settings.translate.nextflow.SUBWORKFLOW_OUTDIR,
        settings.translate.nextflow.CODE_FILES_OUTDIR,
    ]

    file_register: NFFileRegister = NFFileRegister()
    item_register: NFItemRegister = NFItemRegister()

    def __init__(self):
        super().__init__(name="nextflow")

    @classmethod
    def translate_workflow_internal(cls, wf: WorkflowBase) -> Tuple[Any, dict[str, Any]]:

        # set class variables to avoid passing junk params
        settings.translate.nextflow.BASE_OUTDIR = cls.basedir

        # blank scope - main wf has not parent
        scope = Scope()

        # register params and channels for workflow inputs
        params.clear()
        channels.clear()
        # preprocessing.ensure_unique_task_names(wf)
        preprocessing.register_minimal_task_inputs(wf)
        preprocessing.register_params_channels(wf)
        preprocessing.register_data_sources(wf)

        # main logic
        cls.update_files('', '', scope, wf, sources={})

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
        name: str,
        alias: str,
        scope: Scope,
        tool: Workflow | CommandTool | PythonTool,
        sources: dict[str, Any],
        scatter: Optional[ScatterDescription]=None,
        ) -> None:
        
        """
        1. generate all nextflow items (currently: nfgen.NFProcess, nfgen.Workflow, nfgen.ChannelOperation).
        2. write any nextflow items that are processes or workflows to file.
        3. return the nextflow items (so a workflow in scope above can generate process / workflow calls).
        """
        subtype: str = naming.constructs.get_construct_name(tool, scope)

        # groovy code for plumbing
        # TODO: reimplement
        # if scatter and scatter.method == ScatterMethod.cross:
        #     operation_item = nfgen.channels.gen_scatter_cross_operation(sources, scatter)
        #     cls.item_register.add(scope, operation_item)

        # command tool
        if isinstance(tool, CommandTool):

            # groovy library imports & groovy functions used in process
            # item: imports
            # item: functions
            imports_item = parsing.process.gen_imports_for_process(tool)
            functions_item = parsing.process.gen_functions_for_process(tool)
            if imports_item:
                cls.item_register.add(scope, imports_item)
            if functions_item:
                cls.item_register.add(scope, functions_item)

            # item: process
            process_item = parsing.process.gen_process_from_cmdtool(name, alias, tool, sources, scope)
            process_item = cls.handle_container(scope, tool, process_item)
            cls.item_register.add(scope, process_item)
            
            # file: process file
            process_file = NFFile(
                name=name,
                subtype=subtype,
                items=cls.item_register.get(scope),
            )
            cls.file_register.add(scope, process_file)

        # python tool
        elif isinstance(tool, PythonTool):

            # # groovy functions used in process
            # functions_item = nfgen.parsing.process.gen_functions_for_process(tool)
            # if functions_item:
            #     cls.item_register.add(scope, functions_item)
            
            # process
            process_item = parsing.process.gen_process_from_codetool(name, alias, tool, sources, scope)
            process_item = cls.handle_container(scope, tool, process_item)
            cls.item_register.add(scope, process_item)

            # file            
            process_file = NFFile(
                name=name,
                subtype=subtype,
                items=cls.item_register.get(scope),
            )
            cls.file_register.add(scope, process_file)

        # workflow
        elif isinstance(tool, WorkflowBase):  # type: ignore
            
            # parse, update items & files for sub elements (tool / workflow calls)
            wf = tool
            for substep in wf.step_nodes.values():
                scope_copy = deepcopy(scope)
                scope_copy.update(substep)
                cls.update_files(
                    name=naming.constructs.gen_varname_workflow(substep.tool.id()),
                    alias=naming.constructs.gen_varname_workflow(substep.id()),
                    scope=scope_copy,
                    tool=substep.tool,
                    sources=substep.sources,
                    scatter=substep.scatter
                )

            # --- item: imports ---
            # (each file created during parsing sub elements needs to be imported)
            imports_item: Optional[NFImportsBlock] = None
            imports_list: list[NFImport] = []
            for nf_file in cls.file_register.get_children(scope):
                # get the relative import path. this is ugly last-min code. 
                # (main wf)
                if len(scope.labels) == 1: 
                    source = f'./{nf_file.path}'
                # (subwf)
                elif settings.translate.nextflow.SUBWORKFLOW_OUTDIR in nf_file.path:
                    folder = os.path.split(nf_file.path)[-1]
                    source = f'./{folder}'
                elif settings.translate.nextflow.PROCESS_OUTDIR in nf_file.path:
                    source = f'../{nf_file.path}'
                else:
                    raise NotImplementedError
                
                # get the task definition we want to import.
                # only 1 task per file.
                task = [x for x in nf_file.items if isinstance(x, NFProcess) or isinstance(x, NFWorkflow)][0]
                nf_item = NFImportItem(name=task.name, alias=task.alias)
                nf_import = NFImport(method='include', items=[nf_item], source=source)
                imports_list.append(nf_import)

            if imports_list:
                # methods = ['include'] * len(imports_list)
                imports = [i.get_string() for i in imports_list]
                declarations = []  # declarations are for instantiating etc def jsonSlurper = new JsonSlurper()  
                imports_item = NFImportsBlock(imports, declarations)
                cls.item_register.add(scope, imports_item)

            # --- item: functions ---
            # (nothing for now)
            functions_item: Optional[NFFunctionsBlock] = None
            if functions_item:
                cls.item_register.add(scope, functions_item)
            
            # --- item: channels (if in main scope) ---
            channels_item: Optional[NFChannelsBlock] = None
            if len(scope.labels) == 1:
                relevant_channels = channels.getall(scope)
                if relevant_channels:
                    channels_item = NFChannelsBlock(relevant_channels)

            if channels_item:
                cls.item_register.add(scope, channels_item)
 
            # --- item: workflow body ---
            workflow_item = parsing.workflow.gen_workflow(
                scope=scope,
                name=name,
                alias=alias,
                wf=wf,
                item_register=cls.item_register
            )
            cls.item_register.add(scope, workflow_item)

            # ---file: workflow file ---
            # finally, create file for this workflow
            workflow_file = NFFile( 
                name=name,
                subtype=subtype,
                items=cls.item_register.get(scope)
            )
            cls.file_register.add(scope, workflow_file)

        else:
            raise Exception(
                f"Nextflow translation for this {tool.__class__} is not yet supported"
            )

    @classmethod
    def translate_tool_internal(cls, tool: CommandTool) -> str:
        """
        Generate Nextflow process for Janis CommandTool

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        settings.translate.nextflow.MODE = 'tool'
        preprocessing.register_data_sources(tool)
        
        sources: dict[str, Any] = {}
        scope: Scope = Scope()
        
        # hack workaround
        # scope.items = [ToolScopeItem(tool.id())]

        # groovy library imports & groovy functions used in process
        imports_item = parsing.process.gen_imports_for_process(tool)
        functions_item = parsing.process.gen_functions_for_process(tool)
        if imports_item:
            cls.item_register.add(scope, imports_item)
        if functions_item:
            cls.item_register.add(scope, functions_item)

        # process
        name = naming.constructs.gen_varname_process(tool.id())
        alias = None
        process_item = parsing.process.gen_process_from_cmdtool(name, alias, tool, sources, scope)
        process_item = cls.handle_container(scope, tool, process_item)
        cls.item_register.add(scope, process_item)
        
        # file
        process_file = NFFile(
            name=tool.id(),  # TODO here name clash checking
            subtype='tool',
            items=cls.item_register.get(scope),
        )
        settings.translate.nextflow.MODE = 'workflow'
        return process_file.get_string()

    @classmethod
    def translate_code_tool_internal(cls, tool: CodeTool) -> str:
        """
        Generate Nextflow process for Janis CodeTool

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        name = naming.constructs.gen_varname_process(tool.id())
        alias = None
        raise NotImplementedError
        settings.translate.nextflow.MODE = 'tool'
        
        if isinstance(tool, PythonTool):
            preprocessing.register_data_sources(tool)
            
            process = cls.gen_process_from_codetool(tool=tool, sources={}, scope=[])
            process = cls.handle_container(scope, tool, process)

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
            nf_file = NFFile(subtype='', imports=imports, items=items)

            return nf_file.get_string()
        else:
            raise Exception("Only PythonTool code tool is supported for the moment.")
        
    @classmethod
    def handle_container(
        cls,
        scope: Scope,
        tool: Tool,
        process: NFProcess,
    ) -> NFProcess:
        """
        Add container information to a Nexflow NFProcess object

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
        container = None

        if settings.translate.WITH_CONTAINER:
            container = (
                NextflowTranslator.get_container_override_for_tool(tool)
                or tool.container()
            )

        if container is not None:
            container_expr = unwrap_expression(
                val=container, 
                scope=scope, 
            )
            directive = NFContainerDirective(container_expr)
            process.directives.append(directive)

        elif not settings.translate.ALLOW_EMPTY_CONTAINER:
            raise Exception(
                f"The tool '{tool.id()}' did not have a container and no container override was specified. "
                f"Although not recommended, Janis can export empty docker containers with the parameter "
                f"'allow_empty_container=True' or --allow-empty-container"
            )

        return process

    @classmethod
    def translate_helper_files(cls, tool: Tool) -> dict[str, str]:
        """
        Generate a dictionary of helper files to run Nextflow.
        Key of the dictionary is the filename, the value is the file content

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        code_files = cls.gen_python_code_files(tool)
        template_files = cls.gen_template_files(tool)
        helpers = template_files | code_files
        return helpers
    
    @classmethod
    def gen_python_code_files(cls, tool: Tool) -> dict[str, str]:
        # Python files for Python code tools
        files: dict[str, str] = {}

        if isinstance(tool, PythonTool):
            # helpers["__init__.py"] = ""
            #helpers[f"{tool.versioned_id()}.py"] = cls.gen_python_script(tool)
            subdir = settings.translate.nextflow.CODE_FILES_OUTDIR
            filename = f'{tool.id()}.py'
            filepath = os.path.join(subdir, filename)
            files[filepath] = tool.prepared_script(SupportedTranslation.NEXTFLOW)

        elif isinstance(tool, WorkflowBase):
            for step in tool.step_nodes.values():
                step_code_files = cls.gen_python_code_files(step.tool)
                files = files | step_code_files # merge dicts
        
        return files

    @classmethod
    def gen_template_files(cls, tool: Tool) -> dict[str, str]:
        # files from tool.files_to_create
        files: dict[str, str] = {}

        if isinstance(tool, CommandTool):
            if tool.files_to_create():
                for name, contents in tool.files_to_create().items():
                    if not isinstance(name, str):
                        # If name is a File or Directory, the entryname field overrides the value of basename of the File or Directory object 
                        print()
                    
                    if isinstance(contents, str):
                        assert(not name.startswith('unnamed_'))
                        if '<js>' in contents:
                            # ignore, print error message for user
                            pass
                        else:
                            # create file
                            files[name] = contents
                    
                    elif isinstance(contents, InputSelector):
                        tinput_name = contents.input_to_select
                        tinput = tool.inputs_map()[tinput_name]
                        if isinstance(tinput.intype, File | Directory):
                            print('ignored staging File into process')
                        else:
                            print('ignored staging String into process')
                        # # js evaluates to a file: add referenced file to output directory
                        # if name.startswith('unnamed_'):
                        #     # dont override filename
                        #     pass
                        # else:
                        #     # override filename
                        #     pass
                        # print()
                    
                    else:
                        raise NotImplementedError
        
        elif isinstance(tool, WorkflowBase):
            for step in tool.step_nodes.values():
                step_template_files = cls.gen_template_files(step.tool)
                files = files | step_template_files # merge dicts

        return files

    @classmethod
    def unwrap_expression(cls, expression: Any) -> Any:
        pass

    @classmethod
    def build_inputs_dict(
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
        # TODO
        return {}
        raise NotImplementedError
        scope: Scope = Scope()
        if additional_inputs:
            params.register_params_for_additional_inputs(additional_inputs, scope)
        if merge_resources:
            resources_input = cls.build_resources_input(
                tool,
                hints,
                max_cores=max_cores,
                max_mem=max_mem,
                max_duration=max_duration,
            )
            params.register_params_for_resources(resources_input)
        return params.serialize()

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

        build_inputs_dict() and stringify_translated_inputs() must be 
        implemented by any subclass of TranslationBase, but these don't
        properly capture what we want to do for  

        The fed inputs are ignored because we do not want a dict. 
        Rather, we want the registered nfgen.params. We can use their 
        properties to create a nicer config file structure.
        
        :param inputs:
        :type inputs:
        :return:
        :rtype:
        """
        return parsing.config.generate_config()


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
        return 'config'

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
