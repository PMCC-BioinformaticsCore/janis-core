


from abc import ABC, abstractmethod
from enum import Enum
from textwrap import indent
from typing import Optional, Type
from dataclasses import dataclass, field

from janis_core.types import DataType
from janis_core import settings
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType


from .directives import (
    NFProcessDirective,
    NFContainerDirective,
    NFDebugDirective,
    NFCpusDirective,
    NFDiskDirective,
    NFMemoryDirective,
    NFPublishDirDirective,
    NFTimeDirective
)

from .inputs import (
    NFProcessInput, 
    NFPythonToolProcessInput,
    NFTupleProcessInput,
)

from .outputs import NFProcessOutput

INDENT = settings.translate.nextflow.NF_INDENT




class NFProcessScriptType(Enum):
    SCRIPT = "script"
    SHELL = "shell"
    EXEC = "exec"


@dataclass
class NFProcess:
    name: str
    script: str
    script_type: NFProcessScriptType = NFProcessScriptType.SCRIPT
    script_quote: Optional[str] = '"'
    directives: list[NFProcessDirective] = field(default_factory=list)
    inputs: list[NFProcessInput] = field(default_factory=list)
    outputs: list[NFProcessOutput] = field(default_factory=list)
    
    main_exec: Optional[str] = None
    when: Optional[str] = None  # TODO unimplemented?
    pre_script: Optional[str] = None
    
    @property 
    def ordered_directives(self) -> list[NFProcessDirective]:
        return order_nf_directives(self.directives)
    
    @property 
    def ordered_inputs(self) -> list[NFProcessInput]:
        return order_nf_process_inputs(self.inputs)
    
    @property 
    def ordered_outputs(self) -> list[NFProcessOutput]:
        return order_nf_process_outputs(self.outputs)

    @property
    def formatted_directives(self):
        if not self.directives:
            return None
        return "\n".join(INDENT + d.get_string() for d in self.ordered_directives)
    
    @property
    def formatted_inputs(self):
        if not self.inputs:
            return None
        return indent(
            "input:\n" + "\n".join(i.get_string() for i in self.ordered_inputs), 
            INDENT
        )

    @property
    def formatted_outputs(self):
        if not self.outputs:
            return None
        return indent(
            "output:\n" + "\n".join(o.get_string() for o in self.ordered_outputs),
            INDENT,
        )
    
    @property
    def formatted_exec(self) -> Optional[str]:
        if self.main_exec is not None:
            outstr = ''
            outstr += 'exec:'
            if self.main_exec != '':
                outstr += f'\n{self.main_exec}'
            outstr = indent(outstr, INDENT)
            return outstr
        return None
    
    @property
    def formatted_script(self):
        script_body = str(self.script).strip()
        
        if self.script_type == NFProcessScriptType.SCRIPT:
            script = ''
            script += f'{self.script_type.value}:\n'
            script += f'{self.pre_script}\n' if self.pre_script else ''
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script += f'{script_body}\n'
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script = indent(script, INDENT)
        else:
            script = indent(script_body, INDENT)
        
        return script

    def get_string(self) -> str:
        possible_items = [
            self.formatted_directives,
            self.formatted_inputs,
            self.formatted_outputs,
            self.formatted_exec,
            self.formatted_script,
        ]
        final_items = [x for x in possible_items if x is not None]
        tool_definition = "\n\n".join(final_items)
        return f"""\
process {self.name} {{
{tool_definition}
}}
"""



### ORDERING ###

# ORDERING: DIRECTIVES

directive_priorities: dict[Type[NFProcessDirective], int] = {
    NFDebugDirective: 0,
    NFContainerDirective: 1,
    NFPublishDirDirective: 2,
    NFCpusDirective: 9,
    NFDiskDirective: 9,
    NFMemoryDirective: 9,
    NFTimeDirective: 9,
}

class DirectiveOrderer(ABC):
    @abstractmethod
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        ...

class PriorityDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        return sorted(directives, key=lambda x: directive_priorities[type(x)])

class AlphabeticalDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        return sorted(directives, key=lambda x: type(x).__name__)

process_directive_strategies = [
    AlphabeticalDirectiveOrderer,
    PriorityDirectiveOrderer,
]

def order_nf_directives(directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
    for orderer in process_directive_strategies:
        directives = orderer().order(directives)
    return directives


# ORDERING: PROCESS INPUTS

class InputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        ...

class AlphabeticalInputStrategy(InputOrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        return sorted(inputs, key=lambda x: self.input_name(x))

    def input_name(self, pinput: NFProcessInput) -> str:
        if isinstance(pinput, NFTupleProcessInput):
            name = pinput.subnames[0]
        else:
            name = pinput.name
        return name

class MandatoryPriorityInputStrategy(InputOrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        out = sorted(inputs, key=lambda x: x.dtype.optional)
        return out

class TypePriorityInputStrategy(InputOrderingStrategy):

    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        out = sorted(inputs, key=lambda x: self.get_priority(x.dtype))
        return out
    
    def get_priority(self, dtype: DataType) -> int:
        dtt = utils.get_dtt(dtype)
        priorities = {
            DTypeType.SECONDARY_ARRAY: 0,
            DTypeType.SECONDARY: 1,
            DTypeType.FILE_PAIR_ARRAY: 2,
            DTypeType.FILE_PAIR: 3,
            DTypeType.FILE_ARRAY: 4,
            DTypeType.FILE: 5,
            DTypeType.FLAG_ARRAY: 6,
            DTypeType.FLAG: 7,
            DTypeType.GENERIC_ARRAY: 8,
            DTypeType.GENERIC: 9,
        }
        return priorities[dtt]

class PythonToolPriorityInputStrategy(InputOrderingStrategy):
    def order(self, inputs: list[NFProcessInput]) -> list[NFProcessInput]:
        out = sorted(inputs, key=lambda x: isinstance(x, NFPythonToolProcessInput), reverse=True)
        return out

process_input_strategies = [
    AlphabeticalInputStrategy,
    MandatoryPriorityInputStrategy,
    TypePriorityInputStrategy,
    PythonToolPriorityInputStrategy,
]

def order_nf_process_inputs(inputs: list[NFProcessInput]) -> list[NFProcessInput]:
    for strategy in process_input_strategies:
        inputs = strategy().order(inputs)
    return inputs



# ORDERING: PROCESS OUTPUTS

class OutputOrderingStrategy(ABC):
    @abstractmethod
    def order(self, outputs: list[NFProcessOutput]) -> list[NFProcessOutput]:
        ...

class AlphabeticalOutputStrategy(OutputOrderingStrategy):
    def order(self, outputs: list[NFProcessOutput]) -> list[NFProcessOutput]:
        return sorted(outputs, key=lambda x: x.name)

process_output_strategies = [
    AlphabeticalOutputStrategy,
]

def order_nf_process_outputs(outputs: list[NFProcessOutput]) -> list[NFProcessOutput]:
    for strategy in process_output_strategies:
        outputs = strategy().order(outputs)
    return outputs

