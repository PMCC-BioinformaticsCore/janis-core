


from typing import Any, Optional
from textwrap import indent

from janis_core import CommandTool, PythonTool, Workflow
from janis_core.workflow.workflow import StepNode
from janis_core.types import DataType, Stdout

from .. import process
from .. import nfgen_utils
from .. import ordering
from .. import settings

from ..unwrap import unwrap_expression
from ..scope import Scope

from . import trace

from .datatype_mismatch import requires_data_operation
from .datatype_mismatch import handle_data_operation

from .edge_cases import satisfies_edge_case
from .edge_cases import handle_edge_case

from .scatter import is_scatter_relationship
from .scatter import handle_scatter_relationship

NF_INDENT = settings.NF_INDENT


# move to unwrap.py

def get_args(step: StepNode, scope: Scope):
    tool: CommandTool | PythonTool | Workflow   = step.tool     
    sources: dict[str, Any]                     = step.sources
    
    call_args: list[str] = []

    # getting the arg value for each required input 
    valid_input_ids = get_input_ids(tool, sources)
    for name in valid_input_ids:
        if name in sources:
            src = sources[name]

            # get basic arg
            arg = unwrap_expression(
                val=src,
                sources=sources,
            )
            if isinstance(arg, list):
                raise NotImplementedError
                call_args += arg

            # plumbing info
            src_scatter: bool = is_src_scatter(src)
            dest_scatter: bool = is_dest_scatter(name, step)
            srctype: DataType = get_src_type(src)
            desttype: DataType = get_dest_type(tool, name)
            
            # handle scatter relationship
            if name == 'three_prime_adapter_read1':
                print()
                
            # if is_scatter_relationship(src_scatter, dest_scatter):
            #     suffix = handle_scatter_relationship(src_scatter, dest_scatter, src_type, dest_type)
            #     arg = f'{arg}{suffix}'
                        
            # handle datatype relationship
            if satisfies_edge_case(src, desttype, tool):
                suffix = handle_edge_case(src, desttype, tool)
                arg = f'{arg}{suffix}'

            elif requires_data_operation(srctype, desttype, src_scatter, dest_scatter):
                suffix = handle_data_operation(srctype, desttype, src_scatter, dest_scatter)
                arg = f'{arg}{suffix}'

            call_args.append(arg)
            
    # add extra arg in case of python tool - the code file.
    # a param with the same name will have already been created. 
    if isinstance(tool, PythonTool):
        scope_joined = scope.to_string(ignore_base_item=True)
        call_args = [f'params.{scope_joined}.code_file'] + call_args

    return call_args


# def get_args_old(step: StepNode, scope: Scope):
#     tool: CommandTool | PythonTool | Workflow   = step.tool     
#     sources: dict[str, Any]                     = step.sources  
#     scatter: Optional[ScatterDescription]       = step.scatter

#     call_args: list[str] = []
#     input_ids = get_input_ids(tool, sources)

#     # getting the arg value for each required input 
#     for name in input_ids:
#         if name in sources:
#             src = sources[name]
#             scatter_target = True if scatter and name in scatter.fields else False
#             scatter_method = scatter.method if scatter else None
#             res = unwrap_expression(
#                 val=src,
#                 sources=sources,
#                 scatter_target=scatter_target,
#                 scatter_method=scatter_method
#             )
#             if isinstance(res, list):
#                 call_args += res
#             else:
#                 call_args.append(res)
            
#     # add extra arg in case of python tool - the code file.
#     # a param with the same name will have already been created. 
#     if isinstance(tool, PythonTool):
#         scope_joined = scope.to_string(ignore_base_item=True)
#         call_args = [f'params.{scope_joined}.code_file'] + call_args

#     return call_args


# helpers:
# identifying which tool inputs we need to provide an arg for 
def get_input_ids(tool: Workflow | CommandTool | PythonTool, sources: dict[str, Any]) -> list[str]:
    # input ids which we need args for (in correct order)
    if isinstance(tool, Workflow):
        return get_input_ids_workflow(tool, sources)
    else:
        return get_input_ids_tool(tool, sources)

def get_input_ids_workflow(tool: Workflow, sources: dict[str, Any]) -> list[str]:
    # sub Workflow - order via workflow inputs
    # (ignore workflow inputs which don't appear in the original janis step call)
    subwf_ids = set(sources.keys())
    subwf_inputs = nfgen_utils.items_with_id(list(tool.input_nodes.values()), subwf_ids)
    subwf_inputs = ordering.order_workflow_inputs(subwf_inputs)
    return [x.id() for x in subwf_inputs]

def get_input_ids_tool(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> list[str]:
    # CommandTool / PythonTool - order via process inputs 
    process_ids = process.inputs.get_process_inputs(sources)
    process_inputs = nfgen_utils.items_with_id(tool.inputs(), process_ids)
    process_inputs = ordering.order_janis_process_inputs(process_inputs)
    return [x.id() for x in process_inputs]


# helpers:
# identifying types for the data source (upstream wf input or step output)
# and the data destination (tool input)
def get_src_type(src: Any) -> Optional[DataType]:
    # the srctype corresponds to either a workflow input, or step output.
    # scattering doesn't matter. 
    dtype = trace.trace_source_datatype(src)
    if isinstance(dtype, Stdout):
        return dtype.subtype
    else:
        return dtype

def get_dest_type(tool: Workflow | CommandTool | PythonTool, name: str) -> DataType:
    # the desttype corresponds to a tool input. 
    # scattering doesn't matter. 
    tinputs = tool.inputs_map()
    tinp = tinputs[name]
    return tinp.intype

def is_src_scatter(src: Any) -> bool:
    return trace.trace_source_scatter(src)

def is_dest_scatter(name: str, step: StepNode) -> bool:
    if step.scatter and name in step.scatter.fields:
        return True
    return False



# formatting process call text
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


