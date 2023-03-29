


from abc import ABC
from dataclasses import dataclass

from janis_core import CommandTool, PythonTool
from janis_core.workflow.workflow import StepNode

from janis_core import settings

"""
TODO explain this
"""

@dataclass
class ScopeItem(ABC):
    label: str
    subtype: str

@dataclass
class ToolScopeItem(ScopeItem):
    label: str
    subtype: str = 'tool'

@dataclass
class WorkflowScopeItem(ScopeItem):
    label: str
    subtype: str = 'workflow'



class Scope:
    def __init__(self, ):
        self.items: list[ScopeItem] = []
        self.items.append(WorkflowScopeItem(settings.translate.nextflow.NF_MAIN_NAME))
    
    # main way to add
    def update(self, step: StepNode) -> None:
        if isinstance(step.tool, (CommandTool, PythonTool)):
            new_item = ToolScopeItem(step.id())
        else:
            new_item = WorkflowScopeItem(step.id())
        self.items.append(new_item)


    @property 
    def current_entity(self) -> str:
        if not self.labels:
            raise RuntimeError
        return self.labels[-1]            

    @property
    def labels(self) -> list[str]:
        return [x.label for x in self.items]
    
    @property
    def subtypes(self) -> list[str]:
        return [x.subtype for x in self.items]

    def to_string(self, ignore_base_item: bool=False) -> str:
        if ignore_base_item and len(self.labels) > 1:
            out = '.'.join(self.labels[1:])
        # elif ignore_base_item:
        #     out = ''
        elif len(self.labels) >= 1:
            out = '.'.join(self.labels)
        else:
            out = ''
        return out

