





from abc import ABC, abstractmethod
from janis_core import ToolInput
from janis_core.workflow.workflow import InputNode



### COMMAND TOOL INPUTS ###

class ToolInputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        """orders input values and returns ordered list"""
        ...

class AlphabeticalToolInputStrategy(ToolInputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.tag)
        return inputs

class MandatoryPriorityToolInputStrategy(ToolInputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.input_type.optional)
        return inputs

class PositionToolInputStrategy(ToolInputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.position if x.position else 0)
        return inputs


tool_input_strategies = [
    AlphabeticalToolInputStrategy(),
    MandatoryPriorityToolInputStrategy(),
    PositionToolInputStrategy(),
]

def order_tool_inputs(inputs: list[ToolInput]) -> list[ToolInput]:
    for strategy in tool_input_strategies:
        inputs = strategy.order(inputs)
    return inputs

def get_tool_input_positions(inputs: list[ToolInput]) -> dict[str, int]:
    inputs = order_tool_inputs(inputs)
    return {inp.id(): i for i, inp in enumerate(inputs)}



### WORKFLOW INPUTS ###

class WFInputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        """orders input values and returns ordered list"""
        ...

class AlphabeticalWFInputStrategy(WFInputOrderingStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        inputs.sort(key=lambda x: x.id())
        return inputs

class MandatoryPriorityWFInputStrategy(WFInputOrderingStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        inputs.sort(key=lambda x: x.datatype.optional)
        return inputs


workflow_input_strategies = [
    AlphabeticalWFInputStrategy(),
    MandatoryPriorityWFInputStrategy(),
]

def order_workflow_inputs(inputs: list[InputNode]) -> list[InputNode]:
    for strategy in workflow_input_strategies:
        inputs = strategy.order(inputs)
    return inputs

def get_workflow_input_positions(inputs: list[InputNode]) -> dict[str, int]:
    inputs = order_workflow_inputs(inputs)
    return {inp.id(): i for i, inp in enumerate(inputs)}
        
        