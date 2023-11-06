

from typing import Any, Optional
from janis_core import ToolInput, ToolArgument, CommandToolBuilder, DataType, Selector, StringFormatter

from ....nfgen_utils import to_groovy
from ....unwrap import unwrap_expression
from ....variables import VariableManager

from .attributes import get_attributes


def prefix_str(tinput: ToolInput | ToolArgument) -> str:
    if tinput.prefix is not None:
        return tinput.prefix
    return ''

def spacer_str(tinput: ToolInput | ToolArgument) -> str:
    if tinput.separate_value_from_prefix == False:  # type: ignore
        return ''
    else:
        return ' '
    
def delim_str(tinput: ToolInput) -> str:
    if tinput.separator is not None:
        return tinput.separator
    return ' '

def unwrap(val: Any, tool: CommandToolBuilder, vmanager: VariableManager, quote_strings: Optional[bool]=None) -> Any:
    apply_braces = True if isinstance(val, Selector) and not isinstance(val, StringFormatter) else False
    return unwrap_expression(
        val=val,
        context='process_script',
        variable_manager=vmanager,
        tool=tool,
        apply_braces=apply_braces,
        strquote_override=False
    )

def eval_cmdline_targ(
    arg: ToolArgument, 
    tool: CommandToolBuilder, 
    vmanager: VariableManager, 
    shell_quote: Optional[bool]=None
    ) -> str:

    prefix = prefix_str(arg)
    spacer = spacer_str(arg)

    if isinstance(arg.value, list):
        values = [unwrap(elem, tool, vmanager) for elem in arg.value]
        value = spacer.join(values)
    
    else:
        value = unwrap(arg.value, tool, vmanager)

    if shell_quote:
        value = f'"{value}"'
    
    if arg.prefix is not None:
        cmdline = f'{prefix}{spacer}{value}'
    else:
        cmdline = f'{value}'

    if tool.id() == 'BwaMemSamtoolsView' and arg.prefix == '-R':
        print(cmdline)
        print()

    return cmdline


def eval_cmdline_tinput(
    default: Any, 
    tinput: ToolInput, 
    tool: CommandToolBuilder, 
    vmanager: VariableManager,
    apply_prefix: bool,
    quote_strings: Optional[bool]=None
    ) -> str:

    attributes = get_attributes(tinput)

    dtype: DataType = tinput.input_type  # type: ignore
    prefix = prefix_str(tinput)
    spacer = spacer_str(tinput)
    delim = delim_str(tinput)

    if isinstance(default, list) and attributes.prefixeach:
        values: list[str] = []
        for elem in default:
            elem = unwrap(elem, tool, vmanager, quote_strings)
            elem = to_groovy(elem, dtype=dtype)
            elem = f'{prefix}{spacer}{elem}'
            values.append(elem)
        cmdline = f'"{delim.join(values)}"'

    elif isinstance(default, list) and attributes.prefix and apply_prefix:
        values = [unwrap(elem, tool, vmanager, quote_strings) for elem in default]
        values = [to_groovy(elem, dtype=dtype) for elem in values]
        value = delim.join(values)
        cmdline = f'"{prefix}{spacer}{value}"'
    
    elif isinstance(default, list):
        values = [unwrap(elem, tool, vmanager, quote_strings) for elem in default]
        values = [to_groovy(elem, dtype=dtype) for elem in values]
        value = delim.join(values)
        cmdline = f'"{value}"'
    
    elif attributes.prefix and apply_prefix:
        value = unwrap(default, tool, vmanager, quote_strings)
        value = to_groovy(value, dtype=dtype)
        cmdline = f'{prefix}{spacer}{value}'
    
    else:
        value = unwrap(default, tool, vmanager, quote_strings)
        value = to_groovy(value, dtype=dtype)
        cmdline = f'{value}'
    
    return cmdline
