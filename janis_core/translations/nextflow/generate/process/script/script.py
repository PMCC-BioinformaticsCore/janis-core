

from abc import ABC, abstractmethod
from typing import Optional

from janis_core import ToolInput, ToolArgument, CommandTool
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ....variables import VariableManager
from ....variables import VariableHistory
from ....variables import VariableType
from ....variables import Variable

from .ctype import CType, get_ctype
from . import common
from . import selection
from . import ordering
from . import attributes

from janis_core import settings

def gen_script_lines(tool: CommandTool, vmanager: VariableManager) -> list[str]:
    cmdlines: list[str] = []
    # cmdlines += gen_preprocessing_lines(tool, vmanager)
    cmdlines += gen_tool_execution_lines(tool, vmanager)
    return cmdlines

# def gen_preprocessing_lines(tool: CommandTool, vmanager: VariableManager) -> list[str]:
#     cmdlines: list[str] = []
#     pp_inputs = selection.preprocessing_inputs(tool, vmanager)
#     pp_inputs = ordering.script_inputs_arguments(pp_inputs)

#     for tinput in pp_inputs:
#         dtt = utils.get_dtt(tinput.input_type)
        
#         if dtt == DTypeType.SECONDARY_ARRAY:
#             cmdlines += gen_secondary_array_preprocessing(tinput, vmanager)
#         elif dtt == DTypeType.SECONDARY:
#             cmdlines += gen_secondary_preprocessing(tinput, vmanager)
#         else:
#             raise RuntimeError
#     return cmdlines


### BASH PREPROCESSING ###
INDENT = settings.translate.nextflow.NF_INDENT 

# def gen_secondary_array_preprocessing(tinput: ToolInput, vmanager: VariableManager) -> list[str]:
#     header = f'###  ###'

# def gen_secondary_preprocessing(tinput: ToolInput, vmanager: VariableManager) -> list[str]:
#     names = vmanager.get(tinput.id()).original.value
#     primary = names[0]
#     pp: list[str] = []
#     pp.append(f'### preprocessing: grouping secondary files for optional input (workaround) ###\n')
#     pp.append(f'\n')
#     pp.append(f'if [ ${{bam.last()}} != NULL ]\n')
#     pp.append(f'then\n')
#     for name in names:
#         pp.append(f'{INDENT}ln -s ${{{name}}} {primary}/${{{name}.last()}}\n')
#     pp.append('fi\n')
#     pp.append('\n')
#     return pp 

### TOOL EXECUTION ###
    
def gen_tool_execution_lines(tool: CommandTool, vmanager: VariableManager) -> list[str]:
    cmdlines: list[str] = []
    
    all_ins_args = selection.all_script_inputs_arguments(tool, vmanager)
    ordered_ins_args = ordering.script_inputs_arguments(all_ins_args)

    for item in ordered_ins_args:
        
        if selection.is_script_reference_input(item, vmanager):
            formatter = InputReferenceFormatter(item, tool, vmanager) # type: ignore
        elif selection.is_script_value_input(item, vmanager):
            formatter = InputValueFormatter(item, tool, vmanager)     # type: ignore
        elif selection.is_script_argument(item):
            formatter = ArgumentFormatter(item, tool, vmanager)
        else:
            raise RuntimeError
            
        cmdlines += formatter.format()
    
    return cmdlines


### CMDLINE FORMATTING CLASSES ###

class ScriptFormatter(ABC):

    @abstractmethod
    def format(self) -> list[str]:
        ...

### CMDLINE INPUT REFERENCE ###

class InputReferenceFormatter:
    SCRIPT_FMT1 = '${{{src}}}' 
    SCRIPT_FMT2 = '{prefix}{spacer}${{{src}}}'

    def __init__(self, tinput: ToolInput, tool: CommandTool, vmanager: VariableManager) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        self.attributes = attributes.get_attributes(tinput)
        self.ctype = get_ctype(tinput)
        self.dtt = utils.get_dtt(tinput.input_type)

        if self.tool.id() == 'cutadapt' and self.tinput.id() == 'outputPrefix':
            print()

    def format(self) -> list[str]:
        """
        generates a script line for this ToolInput. 
        format is decided by self.itype.
        """
        if self.dtt == DTypeType.FILE_PAIR:
            cmdline = self.SCRIPT_FMT1.format(src=self.cvar)

        elif self.ctype in [
            CType.OPT_BASIC, 
            CType.OPT_BASIC_ARR, 
            CType.OPT_DEFAULT, 
            CType.OPT_DEFAULT_ARR,
        ]:
            cmdline = self.SCRIPT_FMT2.format(
                prefix=common.prefix_str(self.tinput),
                spacer=common.spacer_str(self.tinput),
                src=self.cvar,
            )
        else:
            cmdline = self.SCRIPT_FMT1.format(src=self.cvar)
        
        return [cmdline]
    
    @property
    def varhistory(self) -> VariableHistory:
        return self.vmanager.get(self.tinput.id())
    
    @property
    def cvar(self) -> Variable:
        if self.dtt == DTypeType.FILE_PAIR:
            return self.varhistory.items[1].value
        
        elif self.dtt == DTypeType.SECONDARY:
            if self.attributes.optional:
                return self.varhistory.items[1].value
            else:
                return self.varhistory.items[0].value[0]
        
        else:
            return self.varhistory.current.value



### CMDLINE INPUT VALUE ###
    
class InputValueFormatter(ABC):
    def __init__(self, tinput: ToolInput, tool: CommandTool, vmanager: VariableManager) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        self.attributes = attributes.get_attributes(tinput)
        self.ctype = get_ctype(tinput)
        self.dtt = utils.get_dtt(tinput.input_type)

        if self.dtt in [DTypeType.SECONDARY_ARRAY, DTypeType.SECONDARY]:
            raise RuntimeError

    def format(self) -> list[str]:
        if self.dtt == DTypeType.FLAG:
            cmdline = self.get_cmdline_flag()
        else:
            cmdline = self.get_cmdline_generic()

        lines = [] if cmdline is None else [cmdline]
        return lines

    @property 
    def varhistory(self) -> VariableHistory:
        return self.vmanager.get(self.tinput.id())
        
    def get_cmdline_flag(self) -> Optional[str]:
        if not self.tinput.prefix:
            print()
        assert(self.tinput.prefix)

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
    
    def get_cmdline_generic(self) -> str:
        # tinput not supplied value in any process call, but has default
        if self.varhistory.original.vtype == VariableType.IGNORED and self.tinput.default is not None:
            value = self.tinput.default
        
        # tinput not supplied value in any process call, but is Filename type
        elif self.varhistory.original.vtype == VariableType.IGNORED and utils.is_filename_type(self.tinput.input_type):
            value = self.tinput.input_type
            
        # tinput has consistent static value for each process call
        elif self.varhistory.original.vtype == VariableType.STATIC:
            value = self.varhistory.original.value
        
        else:
            raise RuntimeError

        cmdline = common.eval_cmdline_tinput(
            default=value,
            tinput=self.tinput,
            tool=self.tool,
            vmanager=self.vmanager,
            apply_prefix=True
        )
        return cmdline
    


### CMDLINE ARGUMENT ###

class ArgumentFormatter(ABC):
    def __init__(self, arg: ToolArgument, tool: CommandTool, vmanager: VariableManager) -> None:
        self.tool = tool
        self.arg = arg
        self.vmanager = vmanager

    def format(self) -> list[str]:
        cmdline = common.eval_cmdline_targ(
            arg=self.arg,
            tool=self.tool,
            vmanager=self.vmanager,
            shell_quote=self.arg.shell_quote
        )
        return [cmdline]


    