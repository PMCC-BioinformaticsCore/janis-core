


from typing import Optional

from ...model import XMLParam
from .inputs.Flag import Flag
from .inputs.Option import Option
from .inputs.Positional import Positional
        
from .outputs.RedirectOutput import RedirectOutput
from .outputs.InputOutput import InputOutput
from .outputs.WildcardOutput import WildcardOutput
from .CommandComponent import CommandComponent


def positional(value: Optional[str]=None, gxparam: Optional[XMLParam]=None) -> Positional:
    component = Positional()
    component.gxparam = gxparam
    if value:
        component.values.add(value)
    return component

def flag(prefix: str, gxparam: Optional[XMLParam]=None) -> Flag:
    component = Flag(prefix=prefix)
    component.gxparam = gxparam
    return component

def option( 
    prefix: str,
    gxparam: Optional[XMLParam]=None,
    separator: Optional[str]=None,
    values: Optional[list[str]]=None
) -> Option:
    component = Option(prefix=prefix)
    component.gxparam = gxparam
    if separator:
        component.separator = separator
    if values:
        for val in values:
            component.values.add(val)
    return component

def redirect_output(redirect: str, filepath: str, gxparam: Optional[XMLParam]=None) -> RedirectOutput:
    output = RedirectOutput(redirect, filepath)
    output.gxparam = gxparam
    return output

def input_output(input_component: CommandComponent) -> InputOutput:
    return InputOutput(input_component)

def wildcard_output(gxparam: XMLParam) -> WildcardOutput:
    output = WildcardOutput()
    output.gxparam = gxparam
    return output

def uncertain_output(gxparam: XMLParam) -> WildcardOutput:
    output = WildcardOutput()
    output.gxparam = gxparam
    output.verified = False
    return output
