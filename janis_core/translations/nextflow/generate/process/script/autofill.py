


from typing import Any
from abc import ABC, abstractmethod

from janis_core import ToolInput, CommandTool
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ....variables import VariableHistory
from ....variables import VariableManager
from ....variables import VariableType

from .ctype import CType, get_ctype
from .attributes import get_attributes
from . import common


SCRIPT_FMT1 = '{src}' 
SCRIPT_FMT2 = '{prefix}{spacer}{src}' 



def autofill_script(tinput: ToolInput, tool: CommandTool, vmanager: VariableManager) -> list[str]:
    # when autofill is possible, returns the str expression
    # which can be injected directly into the nf script block.
    # has no pre-script lines. 

    dtt = utils.get_dtt(tinput.input_type)
    autofill_map = {
        DTypeType.SECONDARY_ARRAY: GenericAutoFiller,
        DTypeType.SECONDARY: GenericAutoFiller,
        DTypeType.FILE_PAIR_ARRAY: GenericAutoFiller,
        DTypeType.FILE_PAIR: GenericAutoFiller,
        DTypeType.FILE_ARRAY: GenericAutoFiller,
        DTypeType.FILE: GenericAutoFiller,
        DTypeType.FLAG_ARRAY: FlagArrayAutoFiller,
        DTypeType.FLAG: FlagAutoFiller,
        DTypeType.GENERIC_ARRAY: GenericAutoFiller,
        DTypeType.GENERIC: GenericAutoFiller,
    }

    autofiller = autofill_map[dtt](tool, tinput, vmanager)
    return autofiller.autofill()


class AutoFiller(ABC):
    def __init__(self, tool: CommandTool, tinput: ToolInput, vmanager: VariableManager) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        self.itype = get_ctype(tinput)
        self.dtt = utils.get_dtt(tinput.input_type)
        self.attributes = get_attributes(tinput)
        self.script: list[str] = []

    @abstractmethod
    def autofill(self) -> list[str]:
        ...
    
    @property 
    def varhistory(self) -> VariableHistory:
        return self.vmanager.get(self.tinput.id())
    
    @property
    def value(self) -> Any:
        # ignored input but has default
        if self.varhistory.original.vtype == VariableType.IGNORED:
            value = self.tinput.default
        # static value supplied to input
        elif self.varhistory.original.vtype == VariableType.STATIC:
            value = self.varhistory.original.value
        else:
            raise RuntimeError
        
        return common.eval_default_cmdline(
            default=value,
            tinput=self.tinput,
            tool=self.tool,
            vmanager=self.vmanager,
            apply_prefix=True
        )
    


### ARRAY TYPES ###

class SecondaryArrayAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FilePairArrayAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FileArrayAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FlagArrayAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class GenericArrayAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        cmdline = self.value
        return [cmdline]
    

### SINGLE TYPES ###

class SecondaryAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FilePairAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FileAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        raise NotImplementedError


class FlagAutoFiller(AutoFiller):

    @property
    def value(self) -> Any:
        # ignored input but has default
        if self.varhistory.original.vtype == VariableType.IGNORED:
            value = self.tinput.default
        # static value supplied to input
        elif self.varhistory.original.vtype == VariableType.STATIC:
            value = self.varhistory.original.value
        else:
            raise RuntimeError
        
        value = str(value)
        if value == 'False':
            return None
        else:
            return self.tinput.prefix

    def autofill(self) -> list[str]:
        if self.value is not None:
            cmdline = self.value
            return [cmdline]
        else:
            return []


class GenericAutoFiller(AutoFiller):

    def autofill(self) -> list[str]:
        cmdline = self.value
        return [cmdline]



