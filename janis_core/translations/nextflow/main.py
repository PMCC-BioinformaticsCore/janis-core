from copy import deepcopy
import os
from collections import defaultdict
from typing import Tuple, Optional, Any, Protocol

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

from .model.files.files import NFFile

from .generate.process import generate_processes
from .generate.process import generate_process
from .generate.workflow import generate_workflows
from .generate.files import generate_files 
from .generate.files import generate_file_process 

from . import params
from . import generate
from . import naming
from . import preprocessing

from .scope import Scope




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
    def translate_workflow_internal(cls, wf: Workflow) -> Tuple[Any, dict[str, Any]]:
        # set class variables to avoid passing junk params
        settings.translate.nextflow.BASE_OUTDIR = cls.basedir

        preprocessing.populate_task_inputs_workflowmode(wf, wf)
        processes = generate_processes(wf)
        workflows = generate_workflows(wf, processes)
        files = generate_files(wf, processes, workflows)

        # get the main wf file and all sub files
        main_file = files[wf.id()]
        sub_files = [nf_file for tool_id, nf_file in files.items() if tool_id != wf.id()]

        # return format (gen str for each file)
        main_file_str = main_file.get_string()
        sub_files_str = {sub_file.path: sub_file.get_string() for sub_file in sub_files}
        return (main_file_str, sub_files_str)

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

        preprocessing.populate_task_inputs_toolmode(tool)
        process = generate_process(tool)
        process_file = generate_file_process(process, tool)

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
        assert(isinstance(tool, PythonTool))
        settings.translate.nextflow.MODE = 'tool'

        preprocessing.populate_task_inputs_toolmode(tool)
        process = generate_process(tool)
        process_file = generate_file_process(process, tool)

        return process_file.get_string()

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
        # TODO HERE
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
                        raise NotImplementedError()
                        print()
                    
                    if isinstance(contents, str):
                        assert(not name.startswith('unnamed_'))
                        if '<js>' in contents:
                            # ignore, print error message for user
                            raise NotImplementedError()
                            pass
                        else:
                            # create file
                            path = f'templates/{name}'
                            files[path] = contents
                    
                    elif isinstance(contents, InputSelector):
                        tinput_name = contents.input_to_select
                        tinput = tool.inputs_map()[tinput_name]
                        if isinstance(tinput.intype, File | Directory):
                            raise NotImplementedError()
                            print('ignored staging File into process')
                        else:
                            raise NotImplementedError()
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
        # TODO?
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
        return generate.config.generate_config()

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
