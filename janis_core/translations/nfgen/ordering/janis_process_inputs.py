

from abc import ABC, abstractmethod

from janis_core import ToolInput, TInput
from janis_core.types import File, Array

from .. import nfgen_utils


class OrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
        ...

class AlphabeticalStrategy(OrderingStrategy):
    def order(self, inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
        return sorted(inputs, key=lambda x: x.id())

class PathPriorityStrategy(OrderingStrategy):
    def order(self, inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
        out = sorted(inputs, key=lambda x: self.is_path(x))
        return out
    
    def is_path(self, inp: ToolInput | TInput) -> bool:
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        basetype = nfgen_utils.get_base_type(dtype)
        if isinstance(basetype, File):
            return True
        return False

class TuplePriorityStrategy(OrderingStrategy):
    def order(self, inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
        out = sorted(inputs, key=lambda x: self.is_tuple(x))
        return out
    
    def is_tuple(self, inp: ToolInput | TInput) -> bool:
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        basetype = nfgen_utils.get_base_type(dtype)
        if not isinstance(dtype, Array) and isinstance(basetype, File):
            return True
        return False

class MandatoryPriorityStrategy(OrderingStrategy):
    def order(self, inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
        return sorted(inputs, key=lambda x: self.is_mandatory(x), reverse=True)
    
    def is_mandatory(self, inp: ToolInput | TInput) -> bool:
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        if not dtype.optional:
            return True
        return False


# scatter? 
# array?

process_input_strategies = [
    AlphabeticalStrategy,
    MandatoryPriorityStrategy, 
    PathPriorityStrategy,
    TuplePriorityStrategy,
]

def order_janis_process_inputs(inputs: list[ToolInput | TInput]) -> list[ToolInput | TInput]:
    for strategy in process_input_strategies:
        inputs = strategy().order(inputs)
    return inputs