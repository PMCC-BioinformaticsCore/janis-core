

from typing import Any

from janis_core import Workflow, CommandTool, PythonTool, Tool

from ... import params
from ... import task_inputs

from .populator import TaskInputsPopulatorToolMode
from .populator import TaskInputsPopulatorWorkflowMode
from .common import get_true_workflow_inputs
from ... import naming

# main module funcs

def populate_task_inputs_workflowmode(subwf: Workflow, main_wf: Workflow) -> None:
    """how to populate task inputs when doing workflow translation (workflowmode)"""
    for step in subwf.step_nodes.values():
        # if not already done, formulate task inputs for step task
        if not task_inputs.existsall(step.tool):
            _populate_task_inputs_subtask(step.tool, step.sources, main_wf)
        
        # if subworkflow, do recursively for subworkflow
        if isinstance(step.tool, Workflow):
            populate_task_inputs_workflowmode(step.tool, main_wf)

    # final task is the main workflow
    if subwf.id() == main_wf.id():
        assert(not task_inputs.existsall(main_wf))
        _populate_task_inputs_mainwf(main_wf)

def populate_task_inputs_toolmode(tool: CommandTool | PythonTool) -> None:
    """how to populate task inputs when doing tool translation (toolmode)"""
    populator = TaskInputsPopulatorToolMode(tool)
    populator.populate()


# helper funcs

def _populate_task_inputs_subtask(tool: Tool, sources: dict[str, Any], main_wf: Workflow) -> None:
    """how to populate task inputs for subtask"""
    populator = TaskInputsPopulatorWorkflowMode(tool, sources, main_wf)
    populator.populate()

def _populate_task_inputs_mainwf(wf: Workflow) -> None:
    """
    how to populate task inputs for main wf.
    all the valid workflow inputs start as a param.
    """
    all_tinput_ids = set([x.id() for x in wf.tool_inputs()])
    param_tinput_ids = get_true_workflow_inputs(wf)
    ignored_tinput_ids = all_tinput_ids - param_tinput_ids

    # param inputs
    for tinput_id in param_tinput_ids:
        ti_type = 'param'
        tinput = [x for x in wf.tool_inputs() if x.id() == tinput_id][0]
        subtype = 'main_workflow'
        param = params.register(tinput, task_id=wf.id(), subtype=subtype)
        value = f'params.{param.name}'
        task_inputs.update(wf.id(), ti_type, tinput_id, value)
    
    # ignored inputs
    for tinput_id in ignored_tinput_ids:
        ti_type = 'ignored'
        value = None
        task_inputs.update(wf.id(), ti_type, tinput_id, value)



