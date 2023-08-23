

from janis_core import settings
from janis_core import ToolInput, ToolArgument, CommandTool
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ....variables import VariableManager
from ....variables import VariableType

from .attributes import get_attributes


### PRESCRIPT ###

def prescript_inputs(tool: CommandTool, vmanager: VariableManager) -> list[ToolInput]:
    """
    a ToolInput is necessary when:
    - its VariableType is in [VariableType.TASK_INPUT, VariableType.PARAM, VariableType.STATIC]
    - its VariableType is VariableType.IGNORED and its default value is not None
    - its VariableType is VariableType.IGNORED and it is a Filename type
    """
    out: list[ToolInput] = []
    
    for inp in tool.inputs():
        attributes = get_attributes(inp)
        vtype = vmanager.get(inp.id()).original.vtype
        
        if vtype in [VariableType.TASK_INPUT, VariableType.PARAM]:
            dtt = utils.get_dtt(inp.input_type)  # type: ignore
 
            # these types need special prescript processing
            if dtt in [
                DTypeType.SECONDARY_ARRAY,
                DTypeType.SECONDARY,
                DTypeType.FILE_PAIR_ARRAY,
                DTypeType.FILE_PAIR,
                DTypeType.FILENAME,
                DTypeType.FLAG,
            ]:
                out.append(inp)
            
            # these types always need an array join (which may be inbedded in a conditional statement)
            elif dtt in [
                DTypeType.FILE_ARRAY,
                DTypeType.FLAG_ARRAY,
                DTypeType.GENERIC_ARRAY,
            ]:
                out.append(inp)
            
            # elif these inputs need a conditional statement
            elif attributes.optional or attributes.default:
                out.append(inp)

    return out


### SCRIPT ###

# def preprocessing_inputs(tool: CommandTool, vmanager: VariableManager) -> list[ToolInput]:
#     valid: list[ToolInput] = []
    
#     script_ref_inputs = script_reference_inputs(tool, vmanager)
#     for tinput in script_ref_inputs:
#         dtt = utils.get_dtt(tinput.input_type)
#         attributes = get_attributes(tinput)
        
#         if dtt == DTypeType.SECONDARY_ARRAY and attributes.optional:
#             valid.append(tinput)
#         elif dtt == DTypeType.SECONDARY and attributes.optional:
#             valid.append(tinput)

#     return valid

def all_script_inputs_arguments(tool: CommandTool, vmanager: VariableManager) -> list[ToolInput | ToolArgument]:
    ins_args: list[ToolInput | ToolArgument] = []
    ins_args += script_reference_inputs(tool, vmanager)
    if settings.translate.MODE != 'skeleton':
        ins_args += script_value_inputs(tool, vmanager)
    ins_args += script_arguments(tool)
    return ins_args

def script_reference_inputs(tool: CommandTool, vmanager: VariableManager) -> list[ToolInput]:
    valid: list[ToolInput] = []
    cmdline_tinputs = [x for x in tool.inputs() if x.position is not None or x.prefix is not None]
    for inp in cmdline_tinputs:
        if is_script_reference_input(inp, vmanager):
            valid.append(inp)
    return valid

def is_script_reference_input(inp: ToolInput | ToolArgument, vmanager: VariableManager) -> bool:
    if isinstance(inp, ToolInput):
        vtype = vmanager.get(inp.id()).original.vtype
        if vtype in [VariableType.TASK_INPUT, VariableType.PARAM]:
            return True
    return False

def script_value_inputs(tool: CommandTool, vmanager: VariableManager) -> list[ToolInput]:
    valid: list[ToolInput] = []
    cmdline_tinputs = [x for x in tool.inputs() if x.position is not None or x.prefix is not None]
    for inp in cmdline_tinputs:
        if is_script_value_input(inp, vmanager):
            valid.append(inp)
    return valid

def is_script_value_input(inp: ToolInput | ToolArgument, vmanager: VariableManager) -> bool:
    if isinstance(inp, ToolInput):
        vtype = vmanager.get(inp.id()).original.vtype
        if vtype == VariableType.STATIC:
            return True
        elif vtype == VariableType.IGNORED and inp.default is not None:
            return True
        elif vtype == VariableType.IGNORED and utils.is_filename_type(inp.input_type):
            return True
    return False

def script_arguments(tool: CommandTool) -> list[ToolArgument]:
    if tool.arguments() is not None:
        return [x for x in tool.arguments() if is_script_argument(x)]  # type: ignore
    return []

def is_script_argument(inp: ToolInput | ToolArgument) -> bool:
    # following line seems weird, but ToolInput is actually child class of ToolArgument
    # so have to do it in the negative sense
    if not isinstance(inp, ToolInput):  
        if inp.position is not None or inp.prefix is not None:
            return True
    return False