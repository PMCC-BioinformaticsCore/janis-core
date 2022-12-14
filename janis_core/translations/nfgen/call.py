


from typing import Any, Optional

from janis_core import CommandTool, PythonTool, Workflow
from janis_core.utils.scatter import ScatterDescription

from . import process
from . import nfgen_utils
from . import ordering

from .unwrap import unwrap_expression



def get_args(
    tool: CommandTool | PythonTool | Workflow,
    sources: dict[str, Any],
    scatter: Optional[ScatterDescription],
    ) -> list[str]:
    call_args: list[str] = []

    # subworkflow
    if isinstance(tool, Workflow):
        # sub Workflow - order via workflow inputs
        # (ignore workflow inputs which don't appear in the original janis step call)
        subwf_ids = set(sources.keys())
        subwf_inputs = nfgen_utils.items_with_id(list(tool.input_nodes.values()), subwf_ids)
        subwf_inputs = ordering.order_workflow_inputs(subwf_inputs)
        inputs_ids = [x.id() for x in subwf_inputs]
    
    # everything else
    else:
        # CommandTool / PythonTool - order via process inputs 
        process_ids = process.inputs.get_process_inputs(sources)
        process_inputs = nfgen_utils.items_with_id(tool.inputs(), process_ids)
        process_inputs = ordering.order_janis_process_inputs(process_inputs)
        inputs_ids = [x.id() for x in process_inputs]

    for name in inputs_ids:
        if name in sources:
            src = sources[name]
            scatter_target = True if scatter and name in scatter.fields else False
            scatter_method = scatter.method if scatter else None
            res = unwrap_expression(
                val=src,
                sources=sources,
                scatter_target=scatter_target,
                scatter_method=scatter_method
            )
            if isinstance(res, list):
                call_args += res
            else:
                call_args.append(res)
    return call_args

