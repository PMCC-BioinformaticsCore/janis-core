
from __future__ import annotations

from textwrap import indent
from typing import Optional

from janis_core import settings

from dataclasses import dataclass
from abc import ABC, abstractmethod, abstractproperty

INDENT = settings.translate.nextflow.NF_INDENT



@dataclass
class NFWorkflow(ABC):
    name: str
    alias: Optional[str]
    main: list[str]   

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
workflow {self.name} {{

{self.main_block}

}}
"""


@dataclass
class NFSubWorkflow(NFWorkflow):
    take: list[NFWorkflowTake]
    emit: list[NFWorkflowEmit]

    @property
    def main_block(self) -> str:
        main = "\n".join(self.main)
        main = "main:\n" + main
        return indent(main, INDENT)

    @property
    def take_block(self) -> str:
        return indent(
            "take:\n" + "\n".join(i.get_string() for i in self.take) + '\n', 
            INDENT
        )

    @property
    def emit_block(self) -> str:
        return indent(
            "emit:\n" + "\n".join(i.get_string() for i in self.emit),
            INDENT
        )
    
    def get_string(self) -> str:
        return f"""\
workflow {self.name} {{

{self.take_block}
{self.main_block}
{self.emit_block}

}}
"""


class NFWorkflowTake:
    """
    A nextflow workflow can accept channels as inputs.
    These channels are assigned new names in the workflow `take:` section.
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
    def __init__(self, name: str, as_param: Optional[str] = None):
        self.name = name
        self.as_param = as_param

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