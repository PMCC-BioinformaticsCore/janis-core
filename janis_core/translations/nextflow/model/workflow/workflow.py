
from __future__ import annotations

from textwrap import indent
from typing import Optional

from janis_core import settings
from janis_core.types import DataType, File

from dataclasses import dataclass
from abc import ABC, abstractmethod, abstractproperty
 
INDENT = settings.translate.nextflow.NF_INDENT



@dataclass
class NFWorkflow(ABC):
    name: str
    main: list[str]   
    take: list[NFWorkflowTake]
    emit: list[NFWorkflowEmit]

    @property
    def ordered_take(self) -> list[NFWorkflowTake]:
        return order_workflow_take(self.take)
    
    @property
    def ordered_emit(self) -> list[NFWorkflowEmit]:
        return order_workflow_emit(self.emit)

    @abstractproperty
    def main_block(self) -> str:
        ...

    @abstractmethod
    def get_string(self) -> str:
        ...


@dataclass
class NFMainWorkflow(NFWorkflow):
        
    @property
    def main_block(self) -> str:
        main = "\n".join(self.main)
        return indent(main, INDENT)
    
    def get_string(self) -> str:
        return f"""\
workflow {{

{self.main_block}

}}
"""


@dataclass
class NFSubWorkflow(NFWorkflow):

    @property
    def main_block(self) -> str:
        main = "\n".join(self.main)
        main = "main:\n" + main
        return indent(main, INDENT)

    @property
    def take_block(self) -> str:
        if len(self.take) == 0:
            return ''
        else:
            take = "\ntake:\n" + "\n".join(i.get_string() for i in self.ordered_take) + '\n'
            return indent(take, INDENT)

    @property
    def emit_block(self) -> str:
        if len(self.emit) == 0:
            return ''
        else:
            emit = "emit:\n" + "\n".join(i.get_string() for i in self.ordered_emit)  + '\n\n'
            return indent(emit, INDENT)
    
    def get_string(self) -> str:
        # I hate this
        return f"""\
workflow {self.name} {{
{self.take_block}
{self.main_block}
{self.emit_block}}}
"""


class NFWorkflowTake:
    """
    The only thing binding these channels is their argument position. 
    eg 
        // main wf
        workflow {
            my_pipeline(Channel.of(params.greeting))
        }
        // subwf
        workflow my_pipeline {
            take:
            ch_greeting
        }
    """
    def __init__(self, name: str, tinput_id: str, dtype: DataType):
        self.name = name
        self.tinput_id = tinput_id
        self.dtype = dtype

    def get_string(self) -> str:
        return self.name


class NFWorkflowEmit:
    """
    A nextflow workflow can broadcast channels as outputs. 
    For translation, we always use named outputs.
    eg
        // main wf
        workflow {
            my_pipeline(Channel.of(params.greeting))
            my_pipeline.out.mydata.view()
        }
        // sub wf
        workflow my_pipeline {
            take:
            ch_greeting

            main:
            CONVERTTOUPPER(ch_greeting)

            emit:
            mydata = CONVERTTOUPPER.out.upper
        }
    """
    def __init__(self, name: str, expression: Optional[str] = None):
        self.name = name
        self.expression = expression

    def get_string(self) -> str:
        if self.expression is not None:
            return f"{self.name} = {self.expression}"

        return self.name
    


### ORDERING ###

# ORDERING: TAKE

class TakeSortStrategy(ABC):
    @abstractmethod
    def order(self, take_items: list[NFWorkflowTake]) -> list[NFWorkflowTake]:
        ...

class AlphabeticalTakeSortStrategy(TakeSortStrategy):
    def order(self, take_items: list[NFWorkflowTake]) -> list[NFWorkflowTake]:
        return sorted(take_items, key=lambda x: x.name)

class FileTakeSortStrategy(TakeSortStrategy):
    def order(self, take_items: list[NFWorkflowTake]) -> list[NFWorkflowTake]:
        return sorted(take_items, key=lambda x: isinstance(x.dtype, File), reverse=True)

class MandatoryTakeSortStrategy(TakeSortStrategy):
    def order(self, take_items: list[NFWorkflowTake]) -> list[NFWorkflowTake]:
        return sorted(take_items, key=lambda x: x.dtype.optional == True)

workflow_take_strategies = [
    AlphabeticalTakeSortStrategy, 
    FileTakeSortStrategy,
    MandatoryTakeSortStrategy,
]

def order_workflow_take(take_items: list[NFWorkflowTake]) -> list[NFWorkflowTake]:
    for strategy in workflow_take_strategies:
        take_items = strategy().order(take_items)
    return take_items


# ORDERING: EMIT

class EmitSortStrategy(ABC):
    @abstractmethod
    def order(self, emit_items: list[NFWorkflowEmit]) -> list[NFWorkflowEmit]:
        ...

class AlphabeticalEmitSortStrategy(EmitSortStrategy):
    def order(self, emit_items: list[NFWorkflowEmit]) -> list[NFWorkflowEmit]:
        return sorted(emit_items, key=lambda x: x.name)

workflow_emit_strategies = [
    AlphabeticalEmitSortStrategy, 
]

def order_workflow_emit(emit_items: list[NFWorkflowEmit]) -> list[NFWorkflowEmit]:
    for strategy in workflow_emit_strategies:
        emit_items = strategy().order(emit_items)
    return emit_items
