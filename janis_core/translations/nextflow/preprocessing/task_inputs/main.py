

from typing import Any, Optional

from janis_core import Tool, CommandToolBuilder, PythonTool, WorkflowBuilder
from janis_core import settings

from ... import task_inputs
from .populator import (
    TaskInputsPopulatorToolIngest,
    TaskInputsPopulatorCommandTool,
    TaskInputsPopulatorPythonTool,
    TaskInputsPopulatorSubWorkflow,
    TaskInputsPopulatorMainWorkflow,
)

# main module funcs

def populate_task_inputs(tool: Tool, main_wf: Optional[WorkflowBuilder]=None) -> None:
    assert(isinstance(tool, CommandToolBuilder | PythonTool | WorkflowBuilder))
    if settings.translate.nextflow.ENTITY == 'tool' and isinstance(tool, CommandToolBuilder | PythonTool):
        _populate_task_inputs_toolmode(tool)
    elif settings.translate.nextflow.ENTITY == 'workflow' and isinstance(main_wf, WorkflowBuilder):
        _populate_task_inputs_workflowmode(main_wf, main_wf)
    else:
        raise RuntimeError(f"{tool.id()}: {type(tool)}")

def _populate_task_inputs_workflowmode(subwf: WorkflowBuilder, main_wf: WorkflowBuilder) -> None:
    """how to populate task inputs when doing workflow translation (workflowmode)"""
    for step in subwf.step_nodes.values():
        
        # if not already done, formulate task inputs for step task
        if not task_inputs.existsall(step.tool):
            if isinstance(step.tool, CommandToolBuilder):
                _populate_task_inputs_commandtool(step.tool, step.sources, main_wf)
            elif isinstance(step.tool, PythonTool):
                _populate_task_inputs_pythontool(step.tool, step.sources, main_wf)
            elif isinstance(step.tool, WorkflowBuilder):
                _populate_task_inputs_subwf(step.tool, step.sources, main_wf)
        
        # if subworkflow, do recursively for subworkflow
        if isinstance(step.tool, WorkflowBuilder):
            _populate_task_inputs_workflowmode(step.tool, main_wf)

    # final task is the main workflow
    if subwf.id() == main_wf.id():
        assert(not task_inputs.existsall(main_wf))
        _populate_task_inputs_mainwf(main_wf)


# TOOL MODE 
def _populate_task_inputs_toolmode(tool: CommandToolBuilder | PythonTool) -> None:
    """how to populate task inputs when doing tool translation (toolmode)"""
    populator = TaskInputsPopulatorToolIngest(tool)
    populator.populate()

# WORKFLOW MODE: COMMANDTOOL SUBTASK
def _populate_task_inputs_commandtool(tool: CommandToolBuilder, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
    populator = TaskInputsPopulatorCommandTool(tool, sources, main_wf)
    populator.populate()

# WORKFLOW MODE: PYTHONTOOL SUBTASK
def _populate_task_inputs_pythontool(tool: PythonTool, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
    populator = TaskInputsPopulatorPythonTool(tool, sources, main_wf)
    populator.populate()

# WORKFLOW MODE: WORKFLOW SUBTASK
def _populate_task_inputs_subwf(tool: WorkflowBuilder, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
    populator = TaskInputsPopulatorSubWorkflow(tool, sources, main_wf)
    populator.populate()

# WORKFLOW MODE: MAIN WORKFLOW
def _populate_task_inputs_mainwf(main_wf: WorkflowBuilder) -> None:
    """
    how to populate task inputs for main wf.
    all the valid workflow inputs start as a param.
    """
    populator = TaskInputsPopulatorMainWorkflow(main_wf)
    populator.populate()





