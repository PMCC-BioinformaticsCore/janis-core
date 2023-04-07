


from janis_core import Workflow

from ... import task_inputs

from .populator import TaskInputsPopulator



# main module func

def populate_task_inputs(subwf: Workflow, main_wf: Workflow) -> None:
    for step in subwf.step_nodes.values():
        # if not already done, formulate task inputs for step task
        if not task_inputs.exists(step.tool):
            populator = TaskInputsPopulator(step.tool, step.sources, main_wf)
            populator.populate()
        
        # recursively do for subworkflows 
        if isinstance(step.tool, Workflow):
            populate_task_inputs(step.tool, main_wf)

