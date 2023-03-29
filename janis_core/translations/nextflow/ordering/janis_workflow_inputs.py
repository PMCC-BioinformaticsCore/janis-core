

from abc import ABC, abstractmethod
# from janis_core.workflow.workflow import InputNode
from janis_core.types import File
from janis_core import TInput


class WinpStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[TInput]) -> list[TInput]:
        ...

class AlphabeticalWinpStrategy(WinpStrategy):
    def order(self, inputs: list[TInput]) -> list[TInput]:
        return sorted(inputs, key=lambda x: x.id())

class FileWinpStrategy(WinpStrategy):
    def order(self, inputs: list[TInput]) -> list[TInput]:
        return sorted(inputs, key=lambda x: isinstance(x, File), reverse=True)

class MandatoryWinpStrategy(WinpStrategy):
    def order(self, inputs: list[TInput]) -> list[TInput]:
        return sorted(inputs, key=lambda x: x.intype.optional == True)

workflow_input_strategies = [
    AlphabeticalWinpStrategy, 
    FileWinpStrategy,
    MandatoryWinpStrategy,
]

def order_workflow_inputs(inputs: list[TInput]) -> list[TInput]:
    for strategy in workflow_input_strategies:
        inputs = strategy().order(inputs)
    return inputs