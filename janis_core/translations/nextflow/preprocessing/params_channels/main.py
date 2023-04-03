

from copy import deepcopy
from typing import Any

from janis_core.workflow.workflow import Workflow

from ...scope import Scope
from ... import params

from .ChannelRegistrationHelper import ChannelRegistrationHelper
from .ParamRegistrationHelper import ParamRegistrationHelper
from .helpers_common import get_linkable_params
from .helpers_main import get_param_inputs_to_register_main
from .helpers_main import get_channel_inputs_to_register_main
from .helpers_sub import get_param_inputs_to_register_sub
from .helpers_sub import get_channel_inputs_to_register_sub



def register_params_channels(wf: Workflow) -> None:
    """
    register params & channels for workflow inputs. 
    """
    scope = Scope()
    sources = {}
    do_register_params_channels(wf, sources, scope)

def do_register_params_channels(wf: Workflow, sources: dict[str, Any], scope: Scope):
    
    # link wf InputNodes to params where driven by param in enclosing scope
    # makes params global
    if len(scope.labels) > 1:
        linkable_params = get_linkable_params(wf, sources)
        for uuid, param in linkable_params.items():
            params.create_link(uuid, param)

    param_inputs = get_param_inputs_to_register(wf, sources, scope)
    channel_inputs = get_channel_inputs_to_register(wf, sources, scope)

    for inp in wf.input_nodes.values():
        # params go first
        if inp.id() in param_inputs:
            ParamRegistrationHelper(inp, scope).register()
        if inp.id() in channel_inputs:
            ChannelRegistrationHelper(inp, scope).register()
    
    # repeat for nested workflows (subworkflows)
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        if isinstance(step.tool, Workflow):
            do_register_params_channels(step.tool, step.sources, current_scope)



def get_param_inputs_to_register(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
    if len(scope.labels) > 1:
        return get_param_inputs_to_register_sub(wf, sources, scope)
    else:
        return get_param_inputs_to_register_main(wf, scope)

def get_channel_inputs_to_register(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
    if len(scope.labels) > 1:
        items: set[str] = get_channel_inputs_to_register_sub(wf, sources, scope)
    else:
        items: set[str] = get_channel_inputs_to_register_main(wf, scope)
    return items

