from typing import Any
from janis_core import Workflow

from copy import deepcopy

from ..scope import Scope
from ..plumbing.call import format_process_call
from ..plumbing.call import get_args
from ..process import Process

from .. import channels
from .. import nfgen_utils
from .. import ordering
from .. import settings
from .. import unwrap

from .model import Workflow as NFWorkflow
from .model import WorkflowTake
from .model import WorkflowEmit


def gen_workflow(name: str, scope: Scope, sources: dict[str, Any], wf: Workflow, item_register: Any) -> NFWorkflow:
    """
    Generate a Nextflow Workflow object

    :param workflow:
    :type workflow:
    :param nf_items:
    :type nf_items:
    :param nf_workflow_name:
    :type nf_workflow_name:
    :return:
    :rtype:
    """
    is_subworkflow = True if scope.labels[-1] != settings.NF_MAIN_NAME else False

    take: list[WorkflowTake] = []
    emit: list[WorkflowEmit] = []
    main: list[str] = []

    if is_subworkflow:
        # TAKE
        # which wf inputs should we keep?
        all_inputs = list(wf.input_nodes.values())
        tinput_ids = set(sources.keys())
        tinputs = nfgen_utils.items_with_id(all_inputs, tinput_ids)
        tinputs = ordering.order_workflow_inputs(tinputs)
        
        # confirm channels exist & collect
        chs: list[channels.Channel] = []
        for inp in tinputs:
            assert(channels.exists(inp.uuid))
            ch = channels.get(inp.uuid)
            chs.append(ch)

        # create nf WorkflowTake objects
        for ch in chs:
            take.append(WorkflowTake(ch.name))
        
        # EMIT
        emit: list[WorkflowEmit] = []
        for out in wf.output_nodes.values():
            outname = out.id()
            expression = unwrap.unwrap_expression(val=out.source, in_shell_script=True)
            emit.append(WorkflowEmit(outname, expression))
    
    # MAIN (workflow step calls, channel operations)
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        nf_items = item_register.get(current_scope)
        
        for nf_item in nf_items:
            if isinstance(nf_item, channels.ChannelOperation):
                main.append(nf_item.get_string())
                continue
            elif isinstance(nf_item, Process) or isinstance(nf_item, NFWorkflow):
                entity_name = nf_item.name
            else:
                raise NotImplementedError
        
            args = get_args(step, current_scope)
            main.append(format_process_call(entity_name, args))

    return NFWorkflow(name, main, take, emit, is_subworkflow)