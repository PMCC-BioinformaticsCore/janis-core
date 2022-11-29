


from abc import ABC, abstractmethod
from janis_core import ToolInput
from janis_core.types import File


# essentially the same as above, but has to be different because 
# no shared interface for Workflow and CommandTool (with datatype etc)
class ToolStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        ...

class AlphabeticalToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: x.id())

class FileToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: isinstance(x, File), reverse=True)

class MandatoryToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: x.input_type.optional == True)

tool_input_strategies = [
    #AlphabeticalToolStrategy, 
    FileToolStrategy,
    MandatoryToolStrategy,
]

def order_cmdtool_inputs(inputs: list[ToolInput]) -> list[ToolInput]:
    for strategy in tool_input_strategies:
        inputs = strategy().order(inputs)
    return inputs