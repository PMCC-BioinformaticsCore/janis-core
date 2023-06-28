


from abc import ABC, abstractmethod
from typing import Tuple
from janis_core.ingestion.galaxy.model.workflow.step.inputs import ConnectionInputValue, InputValue, WorkflowInputInputValue

from janis_core.ingestion.galaxy.gx.command.components import Flag
from janis_core.ingestion.galaxy.gx.command.components import Option
from janis_core.ingestion.galaxy.gx.command.components import Positional



def order_positionals(positionals: list[Positional]) -> list[Positional]:
    positionals.sort(key=lambda x: x.cmd_pos)
    return positionals

def order_flags(flags: list[Flag]) -> list[Flag]:
    flags.sort(key=lambda x: x.tag)
    return flags

def order_options(options: list[Option]) -> list[Option]:
    options.sort(key=lambda x: x.tag)
    return options

def order_imports(imports: list[Tuple[str, str]]) -> list[Tuple[str, str]]:
    imports = _order_imports_alphabetical(imports)
    imports = _order_imports_length(imports)
    return imports

def _order_imports_length(imports: list[Tuple[str, str]]) -> list[Tuple[str, str]]:
    imports.sort(key=lambda x: len(x[0] + x[1]), reverse=True)
    return imports

def _order_imports_alphabetical(imports: list[Tuple[str, str]]) -> list[Tuple[str, str]]:
    imports.sort(key=lambda x: f'from {x[0]} import {x[1]}')
    return imports




### INPUT VALUES ###

class InputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        """orders input values and returns ordered list"""
        ...

class AlphabeticalStrategy(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: x.input_tag)
        return invalues

class NotNullPriority(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: x.wrapped_value != 'None', reverse=True)
        return invalues

class PositionalsOptsPositionals(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        opts_pos = self.get_opts_position(invalues)
        top: list[InputValue] = []
        middle: list[InputValue] = []
        bottom: list[InputValue] = []
        for inval in invalues:
            if inval.component and inval.comptype == 'positional':
                if inval.component.cmd_pos < opts_pos:
                    top.append(inval)
                else:
                    bottom.append(inval)
            else:
                middle.append(inval)
        return top + middle + bottom

    def get_opts_position(self, invalues: list[InputValue]) -> int:
        for inval in invalues:
            if inval.component and inval.comptype != 'positional':
                return inval.component.cmd_pos
        return 1

class RuntimeInputPriorityStrategy(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: isinstance(x, WorkflowInputInputValue) and x.is_runtime, reverse=True)
        return invalues

class WorkflowInputPriorityStrategy(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: isinstance(x, WorkflowInputInputValue), reverse=True)
        return invalues

class ConnectionPriorityStrategy(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: isinstance(x, ConnectionInputValue), reverse=True)
        return invalues

class UnlinkedPriorityStrategy(InputOrderingStrategy):
    def order(self, invalues: list[InputValue]) -> list[InputValue]:
        invalues.sort(key=lambda x: x.component is not None)
        return invalues


# changing the order of the objects below changes the 
# ordering priority, as the last ordering method has the highest impact etc

"""
- step tool value order:
    - connection 
    - wflow input 
    - runtime
    - positionals (before opts) 
    - things with values
        - flags 
        - options 
    - things without values
        - flags 
        - options 
    - positionals (after opts)
    - each category above should appear alphabetically
"""

STRATEGIES = [
    AlphabeticalStrategy(),
    NotNullPriority(),
    PositionalsOptsPositionals(),
    RuntimeInputPriorityStrategy(),
    ConnectionPriorityStrategy(),
    WorkflowInputPriorityStrategy(),
    UnlinkedPriorityStrategy()
]

def order_step_inputs(invalues: list[InputValue]):
    for strategy in STRATEGIES:
        invalues = strategy.order(invalues)
    return invalues