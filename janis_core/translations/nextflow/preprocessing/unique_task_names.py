

from copy import deepcopy

from janis_core import Workflow
from janis_core.workflow.workflow import StepNode
from ..scope import Scope


def ensure_unique_task_names(wf: Workflow) -> None:
    scope = Scope()
    task_scopes = gather_task_scopes(wf, scope, {})
    modified_task_scopes = modify_task_scopes(task_scopes)
    rename_step_ids(modified_task_scopes)

def gather_task_scopes(wf: Workflow, scope: Scope, task_scopes: dict[str, StepNode]) -> dict[str, StepNode]:
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        label = current_scope.to_string(ignore_base_item=True)
        task_scopes[label] = step

        if isinstance(step.tool, Workflow):
            task_scopes = gather_task_scopes(step.tool, current_scope, task_scopes)
    
    return task_scopes


def modify_task_scopes(task_scopes: dict[str, StepNode]) -> dict[str, StepNode]:
    modified_task_scopes: dict[str, StepNode] = {}
    
    for label, step in task_scopes.items():
        label_split = label.split('.')
        new_label_split: list[str] = []

        for elem in reversed(label_split):
            new_label_split = [elem] + new_label_split
            new_label = '.'.join(new_label_split)
            if new_label not in modified_task_scopes:
                modified_task_scopes[new_label] = step
                break
    
    return modified_task_scopes


def rename_step_ids(task_scopes: dict[str, StepNode]) -> None:
    for label, step in task_scopes.items():
        label = label.replace('.', '_')
        step.identifier = label