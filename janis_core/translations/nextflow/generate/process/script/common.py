

from typing import Any
from janis_core import ToolInput, CommandTool, DataType

from ....nfgen_utils import to_groovy

from ....unwrap import unwrap_expression
from ....variables import VariableManager

from .attributes import Attributes, get_attributes


def prefix_str(tinput: ToolInput) -> str:
    if tinput.prefix is not None:
        return tinput.prefix
    return ''

def spacer_str(tinput: ToolInput) -> str:
    if tinput.separate_value_from_prefix == False:  # type: ignore
        return ''
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
    vmanager: VariableManager,
    apply_prefix: bool
    ) -> str:

    attributes = get_attributes(tinput)

    dtype: DataType = tinput.input_type  # type: ignore
    prefix = prefix_str(tinput)
    spacer = spacer_str(tinput)
    delim = delim_str(tinput)

    if isinstance(default, list) and attributes.prefixeach:
        values: list[str] = []
        for elem in default:
            elem = unwrap(elem, tool, vmanager)
            elem = to_groovy(elem, dtype=dtype)
            elem = f'{prefix}{spacer}{elem}'
            values.append(elem)
        cmdline = f'"{delim.join(values)}"'

    elif isinstance(default, list) and attributes.prefix and apply_prefix:
        values = [unwrap(elem, tool, vmanager) for elem in default]
        values = [to_groovy(elem, dtype=dtype) for elem in values]
        value = delim.join(values)
        cmdline = f'"{prefix}{spacer}{value}"'
    
    elif isinstance(default, list):
        values = [unwrap(elem, tool, vmanager) for elem in default]
        values = [to_groovy(elem, dtype=dtype) for elem in values]
        value = delim.join(values)
        cmdline = f'"{value}"'
    
    elif attributes.prefix and apply_prefix:
        value = unwrap(default, tool, vmanager)
        value = to_groovy(value, dtype=dtype)
        cmdline = f'{prefix}{spacer}{value}'
    
    else:
        value = unwrap(default, tool, vmanager)
        value = to_groovy(value, dtype=dtype)
        cmdline = f'{value}'
    
    return cmdline
