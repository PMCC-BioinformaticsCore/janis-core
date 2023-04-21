


from typing import Optional

from janis_core import Workflow, PythonTool, CommandTool
from ...model.process import NFProcess
from ...model.workflow import NFWorkflow
from ...model.files import NFFile

from .process import generate_file_process
from .workflow import generate_file_workflow


def generate_files(main_wf: Workflow, nf_processes: dict[str, NFProcess], nf_workflows: dict[str, NFWorkflow]) -> dict[str, NFFile]:
    """generates nextflow files for processes and workflows"""
    nf_files: dict[str, NFFile] = {}
    
    for tool_id, process in nf_processes.items():
        tool = _get_tool(tool_id, main_wf)
        if isinstance(tool, CommandTool | PythonTool):
            nffile = generate_file_process(process, tool)
        else:
            raise Exception(f"Tool '{tool_id}' is not a CommandTool or PythonTool")
        nf_files[tool_id] = nffile
    
    for tool_id, workflow in nf_workflows.items():
        is_subworkflow = True if tool_id != main_wf.id() else False
        tool = _get_tool(tool_id, main_wf)
        assert(isinstance(tool, Workflow))
        nffile = generate_file_workflow(workflow, nf_processes, nf_workflows, tool, is_subworkflow)
        nf_files[tool_id] = nffile

    return nf_files


def _get_tool(tool_id: str, wf: Workflow) -> CommandTool | PythonTool | Workflow:
    """finds & returns workflow tool using tool_id"""
    if wf.id() == tool_id:
        return wf
    
    tool = _do_get_tool(tool_id, wf)
    if not tool:
        raise Exception(f"Tool '{tool_id}' not found in workflow")
    
    return tool

def _do_get_tool(tool_id: str, wf: Workflow) -> Optional[CommandTool | PythonTool | Workflow]:
    """finds & returns workflow tool using tool_id"""
    for step in wf.step_nodes.values():
        if step.tool.id() == tool_id:
            return step.tool

        if isinstance(step.tool, Workflow):
            tool = _do_get_tool(tool_id, step.tool)
            if tool:
                return tool
    
    return None

