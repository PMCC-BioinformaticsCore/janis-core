

from typing import Any
from janis_core import ToolInput, CommandTool, DataType

from ....nfgen_utils import to_groovy

from ....unwrap import unwrap_expression
from ....variables import VariableManager



def prefix_str(tinput: ToolInput) -> str:
    if tinput.prefix is not None:
        return tinput.prefix
    return ''

def spacer_str(tinput: ToolInput) -> str:
    if tinput.separate_value_from_prefix == False:  # type: ignore
        return ''
    elif tinput.separator is not None:
        return tinput.separator
    else:
        return ' '
    
def delim_str(tinput: ToolInput) -> str:
    if tinput.separator is not None:
        return tinput.separator
    return ' '


def unwrap(val: Any, tool: CommandTool, vmanager: VariableManager) -> Any:
    return unwrap_expression(
        val=val,
        context='process_script',
        variable_manager=vmanager,
        tool=tool,
        in_shell_script=True,
    )


def eval_default_cmdline(
    default: Any, 
    tinput: ToolInput, 
    tool: CommandTool, 
    vmanager: VariableManager
    ) -> str:

    dtype: DataType = tinput.input_type  # type: ignore
    prefix = prefix_str(tinput)
    spacer = spacer_str(tinput)
    delim = delim_str(tinput)

    if isinstance(default, list) and tinput.prefix_applies_to_all_elements:
        values: list[str] = []
        for elem in default:
            elem = unwrap(elem, tool, vmanager)
            elem = to_groovy(elem, dtype=dtype)
            elem = f'{prefix}{spacer}{elem}'
            values.append(elem)
        cmdline = ' '.join(values)

    elif isinstance(default, list):
        values = [unwrap(elem, tool, vmanager) for elem in default]
        values = to_groovy(values, dtype=dtype, delim=delim)
        cmdline = f'{prefix}{spacer}{value}'
    
    else:
        value = unwrap(default, tool, vmanager)
        value = to_groovy(value, dtype=dtype)
        cmdline = f'{prefix}{spacer}{value}'
    
    return cmdline
