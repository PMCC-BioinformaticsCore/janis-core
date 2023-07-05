
from abc import ABC, abstractmethod

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from janis_core import ToolInput, ToolArgument
from janis_core.types import Boolean, File, DataType

from .... import nfgen_utils



### PRESCRIPT ORDER ###

class PrescriptOrderingStrategy(ABC):
    @abstractmethod
    def order(self, tinputs: list[ToolInput]) -> list[ToolInput]:
        ...
        
class AlphabeticalPrescriptStrategy(PrescriptOrderingStrategy):
    def order(self, tinputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(tinputs, key=lambda x: x.id())

class DTypePrescriptStrategy(PrescriptOrderingStrategy):
    def order(self, tinputs: list[ToolInput]) -> list[ToolInput]:
        out = sorted(tinputs, key=lambda x: self.get_priority(x.input_type))
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
            DTypeType.FILENAME: 5,
            DTypeType.FLAG_ARRAY: 6,
            DTypeType.FLAG: 7,
            DTypeType.GENERIC_ARRAY: 8,
            DTypeType.GENERIC: 9,
        }
        return priorities[dtt]

prescript_strategies = [
    AlphabeticalPrescriptStrategy,
    DTypePrescriptStrategy,
]

def prescript_inputs(tinputs: list[ToolInput]) -> list[ToolInput]:
    for strategy in prescript_strategies:
        tinputs = strategy().order(tinputs)
    return tinputs




### SCRIPT ORDER ###

class ScriptOrderingStrategy(ABC):
    @abstractmethod
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        ...
        
class AlphabeticalScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: x.prefix or 'zzz')

class ComponentTypeScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        positionals: list[ToolInput | ToolArgument] = []
        flags: list[ToolInput | ToolArgument] = []
        options: list[ToolInput | ToolArgument] = []

        for x in ins_args:
            
            # positionals
            if not x.prefix:
                positionals.append(x)
            
            # flag tool input
            elif isinstance(x, ToolInput) and isinstance(x.input_type, Boolean):
                flags.append(x)

            # opt tool input
            elif isinstance(x, ToolInput):
                options.append(x)
            
            # flag tool arguments
            elif x.value is None:  
                flags.append(x)
            
            # opt tool arguments
            else:
                options.append(x)

        return positionals + options + flags

class InputsPriorityScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: isinstance(x, ToolInput), reverse=True)

class FilePriorityScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        top: list[ToolInput | ToolArgument] = []
        bottom: list[ToolInput | ToolArgument] = []
        for elem in ins_args:
            if isinstance(elem, ToolInput):
                dtype = nfgen_utils.get_base_type_task_input(elem)
                if isinstance(dtype, File):
                    top.append(elem)
                else:
                    bottom.append(elem)
            else:
                bottom.append(elem)
        return top + bottom

class PositionScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda a: (a.position or 0))

class NoPrefixPriorityScriptStrategy(ScriptOrderingStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda a: (a.prefix is None))

script_strategies = [
    # aesthetic
    AlphabeticalScriptStrategy,
    FilePriorityScriptStrategy,
    ComponentTypeScriptStrategy,
    InputsPriorityScriptStrategy,
    # correctness
    NoPrefixPriorityScriptStrategy,
    PositionScriptStrategy
]

def script_inputs_arguments(tinputs: list[ToolInput | ToolArgument]) -> list[ToolInput | ToolArgument]:
    for strategy in script_strategies:
        tinputs = strategy().order(tinputs)
    return tinputs