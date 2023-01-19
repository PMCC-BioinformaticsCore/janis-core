


from typing import Any, Optional
from textwrap import indent

from janis_core import CommandTool, PythonTool, Workflow
from janis_core.workflow.workflow import StepNode
from janis_core.utils.scatter import ScatterDescription

from .. import process
from .. import nfgen_utils
from .. import ordering
from .. import settings

from ..unwrap import unwrap_expression
from ..scope import Scope


NF_INDENT = settings.NF_INDENT

def format_process_call(name: str, inputs: list[str], ind: int=0) -> str:
    if len(inputs) == 0:
        call_str = _call_fmt0(name)
    else:
        call_str = _call_fmt2(name, inputs)
    # elif len(inputs) == 1:
    #     call_str = call_fmt1(name, inputs[0])
    # elif len(inputs) > 1:
        # call_str = call_fmt2(name, inputs)

    return indent(call_str, ind * NF_INDENT)

def _call_fmt0(name: str) -> str:
    return f'{name}()\n'

def _call_fmt1(name: str, input: str) -> str:
    return f'{name}( {input} )\n'

def _call_fmt2(name: str, inputs: list[str]) -> str:
    call_str = f'{name}(\n'
    for i, inp in enumerate(inputs):
        comma = ',' if i < len(inputs) - 1 else ''
        call_str += f'{NF_INDENT}{inp}{comma}\n'
    call_str += ')\n'
    return call_str


def get_args(step: StepNode, scope: Scope):
    tool: CommandTool | PythonTool | Workflow   = step.tool     
    sources: dict[str, Any]                     = step.sources  
    scatter: Optional[ScatterDescription]       = step.scatter

    # input ids which we need args for (in correct order)
    if isinstance(tool, Workflow):
        inputs_ids = get_input_ids_workflow(tool, sources)
    else:
        inputs_ids = get_input_ids_tool(tool, sources)

    call_args: list[str] = []

    # getting the arg value for each required input 
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

    # add extra arg in case of python tool - the code file.
    # a param with the same name will have already been created. 
    if isinstance(tool, PythonTool):
        scope_joined = scope.to_string(ignore_base_item=True)
        call_args = [f'params.{scope_joined}.code_file'] + call_args

    return call_args


def get_input_ids_workflow(tool: Workflow, sources: dict[str, Any]):
    # sub Workflow - order via workflow inputs
    # (ignore workflow inputs which don't appear in the original janis step call)
    subwf_ids = set(sources.keys())
    subwf_inputs = nfgen_utils.items_with_id(list(tool.input_nodes.values()), subwf_ids)
    subwf_inputs = ordering.order_workflow_inputs(subwf_inputs)
    return [x.id() for x in subwf_inputs]


def get_input_ids_tool(tool: CommandTool | PythonTool, sources: dict[str, Any]):
    # CommandTool / PythonTool - order via process inputs 
    process_ids = process.inputs.get_process_inputs(sources)
    process_inputs = nfgen_utils.items_with_id(tool.inputs(), process_ids)
    process_inputs = ordering.order_janis_process_inputs(process_inputs)
    return [x.id() for x in process_inputs]



