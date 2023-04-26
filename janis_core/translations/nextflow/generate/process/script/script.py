

from janis_core import ToolInput, CommandTool

from ....variables import VariableManager
from ....variables import VariableHistory
from ....variables import VariableType

from .ctype import CType, get_ctype
from . import common
from . import autofill
from . import attributes



def gen_script_lines(
    tinput: ToolInput,
    tool: CommandTool,
    vmanager: VariableManager,
) -> list[str]:
    
    if _should_ignore(tinput, vmanager):
        return []
    
    elif _should_autofill(tinput, vmanager):
        return autofill.autofill_script(tinput, tool, vmanager)
    
    else:
        formatter = ScriptFormatter(tinput, tool, vmanager)
        return formatter.format()


### HELPER METHODS ###

def _should_ignore(tinput: ToolInput, vmanager: VariableManager) -> bool:
    # tinput not supplied value in any process call, no default
    # no prescript, no script
    varhistory = vmanager.get(tinput.id())
    if varhistory.original.vtype == VariableType.IGNORED:
        if tinput.default is None:
            return True
    return False

def _should_autofill(tinput: ToolInput, vmanager: VariableManager) -> bool:
    """
    For ToolInputs which are not fed via a process input or a param,
    should a static value be evaluated? 
    """
    varhistory = vmanager.get(tinput.id())

    # tinput not supplied value in any process call, but has default
    # no prescript, no script
    if varhistory.original.vtype == VariableType.IGNORED:
        if tinput.default is not None:
            return True
    
    # tinput has consistent static value for each process call
    # no prescript, script autofilled
    elif varhistory.original.vtype == VariableType.STATIC:
        return True
    return False


### MAIN FORMATTING CLASS ###

SCRIPT_FMT1 = '${{{src}}}' 
SCRIPT_FMT2 = '{prefix}{spacer}${{{src}}}'

class ScriptFormatter:
    def __init__(
        self, 
        tinput: ToolInput,
        tool: CommandTool,
        vmanager: VariableManager,
    ) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        self.attributes = attributes.get_attributes(tinput)
        self.ctype = get_ctype(tinput)

    def format(self) -> list[str]:
        """
        generates a script line for this ToolInput. 
        format is decided by self.itype.
        """
        if self.ctype in [
            CType.OPT_BASIC,
            CType.OPT_BASIC_ARR,
            CType.OPT_DEFAULT,
            CType.OPT_DEFAULT_ARR,
        ]:
            cmdline = SCRIPT_FMT2.format(
                prefix=self.prefix_str,
                spacer=self.spacer_str,
                src=self.varhistory.current.value,
            )
        
        else:
            cmdline = SCRIPT_FMT1.format(
                src=self.varhistory.current.value
            )
        
        return [cmdline]
    
    @property
    def varhistory(self) -> VariableHistory:
        return self.vmanager.get(self.tinput.id())

    @property
    def prefix_str(self) -> str:
        return common.prefix_str(self.tinput)
    
    @property
    def spacer_str(self) -> str:
        return common.spacer_str(self.tinput)
    
    @property
    def delim_str(self) -> str:
        return common.delim_str(self.tinput)



        

    
    