
from __future__ import annotations

from textwrap import indent
from typing import Optional

from janis_core import settings
from janis_core.types import DataType

from dataclasses import dataclass
from abc import ABC, abstractmethod, abstractproperty
 
INDENT = settings.translate.nextflow.NF_INDENT



@dataclass
class NFWorkflow(ABC):
    name: str
    main: list[str]   
    take: list[NFWorkflowTake]
    emit: list[NFWorkflowEmit]

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
            take = "\ntake:\n" + "\n".join(i.get_string() for i in self.take) + '\n'
            return indent(take, INDENT)

    @property
    def emit_block(self) -> str:
        if len(self.emit) == 0:
            return ''
        else:
            emit = "emit:\n" + "\n".join(i.get_string() for i in self.emit)  + '\n\n'
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