





from abc import ABC, abstractmethod
from janis_core import ToolInput


### INPUT VALUES ###

class InputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        """orders input values and returns ordered list"""
        ...

class AlphabeticalStrategy(InputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.tag)
        return inputs

class MandatoryPriorityStrategy(InputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.input_type.optional)
        return inputs

class PositionStrategy(InputOrderingStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        inputs.sort(key=lambda x: x.position if x.position else 0)
        return inputs


strategies = [
    AlphabeticalStrategy(),
    MandatoryPriorityStrategy(),
    PositionStrategy(),
]

def order_tool_inputs(inputs: list[ToolInput]) -> list[ToolInput]:
    for strategy in strategies:
        inputs = strategy.order(inputs)
    return inputs

def get_tool_input_ordering(inputs: list[ToolInput]) -> dict[str, int]:
    inputs = order_tool_inputs(inputs)
    return {inp.tag: i for i, inp in enumerate(inputs)}
        
        