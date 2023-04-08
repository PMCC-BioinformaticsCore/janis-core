

from abc import ABC, abstractmethod

from ..model.process.inputs import (
    NFProcessInput,
    NFPathProcessInput,
    NFTupleProcessInput,
)


class OrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        ...

class AlphabeticalStrategy(OrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        return sorted(inputs, key=lambda x: self.input_name(x))

    def input_name(self, pinput: NFProcessInput) -> str:
        if isinstance(pinput, NFTupleProcessInput):
            name = pinput.subnames[0]
        else:
            name = pinput.name
        return name

class PathPriorityStrategy(OrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        out = sorted(inputs, key=lambda x: isinstance(x, NFPathProcessInput), reverse=True)
        return out

class TuplePriorityStrategy(OrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        out = sorted(inputs, key=lambda x: isinstance(x, NFTupleProcessInput), reverse=True)
        return out


process_input_strategies = [
    AlphabeticalStrategy,
    PathPriorityStrategy,
    TuplePriorityStrategy,
]

def order_nf_process_inputs(inputs: list[NFProcessInput]) -> list[NFProcessInput]:
    for strategy in process_input_strategies:
        inputs = strategy().order(inputs)
    return inputs