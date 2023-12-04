
import os
from typing import Tuple, Any, Optional

from janis_core.tool.commandtool import Tool
from janis_core.code.codetool import CodeTool
from janis_core.code.pythontool import PythonTool
from janis_core.tool.commandtool import CommandToolBuilder
from janis_core.workflow.workflow import Workflow, WorkflowBase, WorkflowBuilder
from janis_core.translations.translationbase import TranslatorBase
from janis_core.translation_deps.supportedtranslations import SupportedTranslation
from janis_core import InputSelector, File, Directory
from janis_core import settings

from .model.workflow import NFWorkflow
from .model.process import NFProcess

from .generate.process import generate_processes
from .generate.process import generate_process
from .generate.workflow import generate_workflows
# from .generate.files import generate_files 
from .generate.files import generate_file_process 
from .generate.files import generate_file_workflow 
from .casefmt import to_case
from . import generate
from . import preprocessing


# HELPERS
def _get_tool(tool_id: str, wf: WorkflowBuilder) -> CommandToolBuilder | CodeTool:
    """finds & returns workflow tool using tool_id"""
    tool = _do_get_tool(tool_id, wf)
    if not tool:
        raise Exception(f"Tool '{tool_id}' not found in workflow")
    return tool

def _do_get_tool(tool_id: str, wf: WorkflowBuilder) -> Optional[CommandToolBuilder | CodeTool]:
    """finds & returns workflow tool using tool_id"""
    for step in wf.step_nodes.values():
        if step.tool.id() == tool_id:
            return step.tool

        if isinstance(step.tool, WorkflowBuilder):
            tool = _do_get_tool(tool_id, step.tool)
            if tool:
                return tool
    return None

def _get_workflow(workflow_id: str, main_wf: WorkflowBuilder) -> WorkflowBuilder:
    if main_wf.id() == workflow_id:
        return main_wf
    workflow = _do_get_workflow(workflow_id, main_wf)
    if not workflow:
        raise Exception(f"Workflow '{workflow_id}' not found in workflow")
    return workflow

def _do_get_workflow(workflow_id: str, wf: WorkflowBuilder) -> Optional[WorkflowBuilder]:
    """finds & returns workflow using workflow_id"""
    for step in wf.step_nodes.values():
        if step.tool.id() == workflow_id:
            return step.tool

        if isinstance(step.tool, Workflow):
            workflow = _do_get_workflow(workflow_id, step.tool)
            if workflow:
                return workflow
    
    return None


class NextflowTranslator(TranslatorBase):
    
    OUTDIR_STRUCTURE: dict[str, str | None] = {
        'main': None,
        'subworkflows': 'subworkflows',
        'tools': 'modules',
        'inputs': None,
        'helpers': 'templates',
        'resources': None,
    }

    def __init__(self):
        super().__init__(name="nextflow")


    ### JANIS -> OUTPUT MODEL MAPPING ###

    def translate_workflow_internal(self, wf: Workflow) -> None:
        # set class variables to avoid passing junk data
        assert(isinstance(wf, WorkflowBuilder))
        settings.translate.nextflow.BASE_OUTDIR = self.basedir

        preprocessing.populate_task_inputs(wf, wf)
        processes = generate_processes(wf)
        workflows = generate_workflows(wf, processes)
        
        for tool_id, nf_process in processes.items():
            tool = _get_tool(tool_id, wf)
            self.add_tool(tool, nf_process)
        
        for wf_id, nf_workflow in workflows.items():
            is_subworkflow = True if wf_id != wf.id() else False
            janis_workflow = _get_workflow(wf_id, wf)
            if is_subworkflow:
                self.add_subworkflow(janis_workflow, nf_workflow)
            else:
                self.main = (janis_workflow, nf_workflow)

    def translate_tool_internal(self, tool: CommandToolBuilder) -> None:
        """
        Generate Nextflow process for Janis CommandToolBuilder

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        assert(isinstance(tool, CommandToolBuilder))
        # check we haven't already translated this
        if self.get_tool(tool) is not None:
            return
        settings.translate.nextflow.ENTITY = 'tool'
        preprocessing.populate_task_inputs(tool)
        tool_nxf = generate_process(tool)
        self.add_tool(tool, tool_nxf)

    def translate_code_tool_internal(self, tool: CodeTool) -> None:
        """
        Generate Nextflow process for Janis CodeTool

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        assert(isinstance(tool, PythonTool))
        # check we haven't already translated this
        if self.get_tool(tool) is not None:
            return
        settings.translate.nextflow.ENTITY = 'tool'
        preprocessing.populate_task_inputs(tool)
        tool_nxf = generate_process(tool)
        self.add_tool(tool, tool_nxf)
    
    def build_inputs_file(self, entity: WorkflowBuilder | CommandToolBuilder | CodeTool) -> None:
        """delegated to external function due to size / complexity"""
        self.inputs_file = generate.config.generate_config()
    
    def build_resources_file(self, entity: WorkflowBuilder | CommandToolBuilder | CodeTool) -> None:
        """
        Always combined for nextflow.
        """
        self.resources_file = None

    @classmethod
    def unwrap_expression(cls, expression: Any) -> Any:
        pass


    ### STRINGIFY ###

    def stringify_translated_workflow(self, internal: WorkflowBuilder, translated: NFWorkflow) -> str:
        # need to get all translated workflows & processes to generate imports. 
        # Can't remember why it needs to be dict - oh well
        
        # main workflow
        nf_workflows = {}
        assert self.main is not None
        intnl_main, trans_main = self.main
        nf_workflows[intnl_main.id()] = trans_main
        
        # subworkflows
        for intnl_subwf, nf_subwf in self.subworkflows:
            nf_workflows[intnl_subwf.id()] = nf_subwf
        
        # processes
        nf_processes = {}
        for intnl_tool, nf_process in self.tools:
            nf_processes[intnl_tool.id()] = nf_process

        is_subworkflow = True if internal.id() == self.main[0].id() else False
        nffile = generate_file_workflow(
            translated, 
            nf_processes, 
            nf_workflows, 
            internal, 
            is_subworkflow
        )
        return nffile.get_string()

    def stringify_translated_tool(self, internal: CommandToolBuilder | CodeTool, translated: NFProcess) -> str:
        nffile = generate_file_process(translated, internal)
        return nffile.get_string()

    
    ### FILENAMES ###

    @staticmethod
    def workflow_filename(workflow: WorkflowBuilder, is_main: Optional[bool]=False) -> str:
        """
        Generate the main workflow filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        if is_main:
            return settings.translate.nextflow.MAIN_WORKFLOW_NAME
        basename = to_case(workflow.id(), settings.translate.nextflow.NF_FILE_CASE)
        return basename + ".nf"

    @staticmethod
    def inputs_filename(workflow: Workflow) -> str:
        """
        Generate the input filename

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        return settings.translate.nextflow.CONFIG_FILENAME

    @staticmethod
    def tool_filename(tool: str | Tool) -> str:
        toolname = tool if isinstance(tool, str) else tool.id()
        basename = to_case(toolname, settings.translate.nextflow.NF_FILE_CASE)
        return basename + ".nf"

    @staticmethod
    def resources_filename(workflow: Workflow) -> str:
        """
        Generate resoureces filename 
        This should never be called.

        :param workflow:
        :type workflow:
        :return:
        :rtype:
        """
        raise RuntimeError
        # return workflow.id() + "-resources.json"

    
    ### VALIDATION ###

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass
