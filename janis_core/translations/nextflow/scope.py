


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
    def __init__(self):
        self.items: list[ScopeItem] = []
        self.items.append(WorkflowScopeItem(settings.translate.nextflow.NF_MAIN_NAME))
    
    def update(self, step: StepNode) -> None:
        new_item = self.create_item(step)
        self.items.append(new_item)

    def create_item(self, step: StepNode) -> ScopeItem:
        if isinstance(step.tool, (CommandTool, PythonTool)):
            return ToolScopeItem(step.id())
        else:
            return WorkflowScopeItem(step.id())

    @property
    def labels(self) -> list[str]:
        return [x.label for x in self.items]
    
    @property
    def subtypes(self) -> list[str]:
        return [x.subtype for x in self.items]

    def to_string(self, ignore_base_item: bool=False) -> str:
        if ignore_base_item and len(self.labels) > 1:
            out = '.'.join(self.labels[1:])
        elif len(self.labels) >= 1:
            out = '.'.join(self.labels)
        else:
            out = ''
        return out

