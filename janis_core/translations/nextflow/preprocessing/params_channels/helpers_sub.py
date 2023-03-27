

from typing import Any

from janis_core.workflow.workflow import Workflow, InputNode, StepNode
from janis_core import translation_utils as utils

from ... import channels
from ... import params
from ...scope import Scope

from .helpers_common import get_true_workflow_inputs



def get_param_inputs_to_register_sub(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
    """
    which sub workflow InputNodes qualify as params we should register?
    in sub workflow, every data source gets a param EXCEPT:
    - if the data source will be a nextflow channel
    - if the data source can be backtraced to a param which already exists. 
    - if the data source is static (ie 1 or 'conservative' etc)
    
    eg we have a global param, which is accessed within a subworkflow.
    """
    all_inputs = get_true_workflow_inputs(wf)
    channel_inputs = get_channel_inputs_to_register_sub(wf, sources, scope)
    global_param_inputs = get_global_param_inputs(wf)
    
    surviving_inputs = all_inputs
    surviving_inputs = surviving_inputs - channel_inputs
    surviving_inputs = surviving_inputs - global_param_inputs
    return surviving_inputs

def get_global_param_inputs(wf: Workflow) -> set[str]:
    all_inputs = get_true_workflow_inputs(wf)
    input_nodes = [wf.input_nodes[x] for x in all_inputs]

    out: set[str] = set()
    for inp in input_nodes:
        if params.exists(inp.uuid):
            out.add(inp.id())
    return out

def get_channel_inputs_to_register_sub(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
    """
    Get the wf inputs for which we will create a nf channel in a sub wf.
    inputs which:
    - are fed from a channel
    - are fed from a step output
    will be channels. others will be params, or have a static value which will be used in the subworkflow. 
    """
    channel_inputs = get_channel_inputs_sub(wf, sources, scope)
    connection_inputs = get_connection_inputs_sub(wf, sources)
    final_inputs = channel_inputs | connection_inputs
    return final_inputs

def get_channel_inputs_sub(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
    out: set[str] = set()
    for name, inp in sources.items():
        node = utils.resolve_node(inp)
        if isinstance(node, InputNode):
            if channels.exists(scope, node.uuid, for_parent=True): 
                out.add(name)
    return out

def get_connection_inputs_sub(wf: Workflow, sources: dict[str, Any]) -> set[str]:
    out: set[str] = set()
    for name, inp in sources.items():
        node = utils.resolve_node(inp)
        if isinstance(node, StepNode):
            out.add(name)
    return out