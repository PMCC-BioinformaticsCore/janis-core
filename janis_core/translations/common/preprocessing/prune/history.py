

from typing import Any
from dataclasses import dataclass, field

from janis_core import WorkflowBuilder, TInput, CommandToolBuilder, PythonTool
from janis_core.workflow.workflow import StepNode
from janis_core.operators.selectors import InputNodeSelector
from janis_core.operators.selectors import StepOutputSelector
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType


@dataclass
class InputTaskInput:
    step_id: str
    value: InputNodeSelector

@dataclass
class ConnectionTaskInput:
    step_id: str
    value: StepOutputSelector

@dataclass
class OtherTaskInput:
    step_id: str
    value: Any

TaskInput = InputTaskInput | ConnectionTaskInput | OtherTaskInput


@dataclass
class TaskInputHistory:
    tinput: TInput
    sources: list[TaskInput] = field(default_factory=list)

    @property
    def is_optional(self) -> bool:
        dtt = utils.get_dtt(self.tinput.intype)
        if dtt == DTypeType.FILENAME:
            return False
        elif self.tinput.intype.optional == True:
            return True
        return False
    
    @property
    def input_sources(self) -> list[InputTaskInput]:
        return [x for x in self.sources if isinstance(x, InputTaskInput)]
    
    @property
    def connection_sources(self) -> list[ConnectionTaskInput]:
        return [x for x in self.sources if isinstance(x, ConnectionTaskInput)]
    
    @property
    def other_sources(self) -> list[OtherTaskInput]:
        return [x for x in self.sources if isinstance(x, OtherTaskInput)]
    
    @property
    def mandatory_input_sources(self) -> list[InputTaskInput]:
        sources = self.input_sources
        sources = [x for x in sources if x.value.input_node.datatype.optional == False]
        return sources
    
    @property
    def placeholder_sources(self) -> list[InputTaskInput]:
        sources: list[InputTaskInput] = []
        for source in self.input_sources:
            # does the source input node source look like a placeholder?
            node = source.value.input_node
            step_id = source.step_id
            if utils.looks_like_placeholder_node(node, step_id, self.tinput.id(), self.tinput.intype):
                sources.append(source)
        return sources
    
    def add_value(self, step_id: str, src: Any) -> None:
        if isinstance(src, InputNodeSelector):
            ti = InputTaskInput(step_id, src)
        elif isinstance(src, StepOutputSelector):
            ti = ConnectionTaskInput(step_id, src)
        else:
            ti = OtherTaskInput(step_id, src)
        self.sources.append(ti)



class TaskInputCollector:
    """
    for a given tool_id, searches the workflow for each step calling that tool.
    records the values provided to each TInput in that step call. 
    """
    def __init__(self, tool: CommandToolBuilder | PythonTool | WorkflowBuilder) -> None:
        self.tool = tool
        self.histories: dict[str, TaskInputHistory] = {}
        self.step_count: int = 0

    def collect(self, wf: WorkflowBuilder) -> None:
        # iterate through workflow steps, finding those which call self.tool
        for step in wf.step_nodes.values():
            if not isinstance(step.tool, CommandToolBuilder | PythonTool | WorkflowBuilder):
                raise RuntimeError
            
            # collect task inputs if step calls self.tool
            if step.tool.id() == self.tool.id():
                inputs_dict: dict[str, Any] = {}
                self.step_count += 1
                
                # update task inputs for tinputs with item in sources
                inputs_dict = inputs_dict | self.gather_sources(step.sources)
                
                # update the dict
                self.update_histories(inputs_dict, step)

            # recursive for nested workflows
            if isinstance(step.tool, WorkflowBuilder):
                self.collect(step.tool)
    
    def gather_sources(self, sources: dict[str, Any]) -> dict[str, InputNodeSelector]:
        out: dict[str, InputNodeSelector] = {}
        for tinput_id, src in sources.items():
            selector = src.source_map[0].source
            out[tinput_id] = selector
        return out

    def update_histories(self, inputs_dict: dict[str, Any], step: StepNode) -> None:
        # update the record of each TInput's history
        for tinput_id, src in inputs_dict.items():

            # add a TaskInputHistory for TInput if not exists
            if tinput_id not in self.histories:
                tinput = [x for x in step.tool.tool_inputs() if x.id() == tinput_id][0]
                history = TaskInputHistory(tinput)
                self.histories[tinput_id] = history
            
            # add a value to this TInput's TaskInputHistory
            if src is not None:
                self.histories[tinput_id].add_value(step.id(), src)

