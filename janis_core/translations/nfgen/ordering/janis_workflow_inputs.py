

from abc import ABC, abstractmethod
from janis_core.workflow.workflow import InputNode
from janis_core.types import File


class WinpStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        ...

class AlphabeticalWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: x.id())

class FileWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: isinstance(x, File), reverse=True)

class MandatoryWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: x.datatype.optional == True)

workflow_input_strategies = [
    AlphabeticalWinpStrategy, 
    FileWinpStrategy,
    MandatoryWinpStrategy,
]

def order_workflow_inputs(inputs: list[InputNode]) -> list[InputNode]:
    for strategy in workflow_input_strategies:
        inputs = strategy().order(inputs)
    return inputs