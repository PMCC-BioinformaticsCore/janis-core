


from typing import Any, Optional

from janis_core import CommandTool, PythonTool, Workflow
from janis_core.utils.scatter import ScatterDescription

from . import process
from . import nfgen_utils
from . import ordering

from .unwrap import unwrap_expression

# from janis_core.graph.steptaginput import StepTagInput

# from janis_core.operators.operator import IndexOperator
# from janis_core.operators.standard import FirstOperator
# from janis_core.operators.selectors import InputNodeSelector, StepOutputSelector, AliasSelector

# from janis_core.workflow.workflow import StepNode, InputNode
# from janis_core.types import DataType
# from janis_core.utils.scatter import ScatterMethod

# from . import settings
# from . import channels
# from . import params
# from .scatter import cartesian_cross_subname
# from .casefmt import to_case


def get_args(
    tool: CommandTool | PythonTool | Workflow,
    sources: dict[str, Any],
    scatter: Optional[ScatterDescription],
    ) -> list[str]:
    call_args: list[str] = []

    # subworkflow
    if isinstance(tool, Workflow):
        process_inputs_names = set(sources.keys())
    
    # everything else
    else:
        process_ids = process.get_process_inputs(sources)
        process_inputs = nfgen_utils.items_with_id(tool.inputs(), process_ids)
        process_inputs = ordering.order_janis_process_inputs(process_inputs)
        process_inputs_names = [x.id() for x in process_inputs]


    for name in process_inputs_names:
        if name in sources:
            src = sources[name]
            scatter_target = True if scatter and name in scatter.fields else False
            scatter_method = scatter.method if scatter else None
            args = unwrap_expression(
                value=src,
                sources=sources,
                scatter_target=scatter_target,
                scatter_method=scatter_method
            )
            call_args += args
    return call_args


"""
unsure if needed?

if settings.PYTHON_CODE_FILE_PATH_PARAM in value:
    path_to_python_code_file = posixpath.join(
        #"$baseDir", cls.DIR_TOOLS, f"{tool.versioned_id()}.py"
        "$baseDir", cls.DIR_TOOLS, f"{tool.id()}.py"
    )
    value = value.replace(
        settings.PYTHON_CODE_FILE_PATH_PARAM, f'"{path_to_python_code_file}"'
    )
"""

# def unwrap(tinput_name: str, node: Any, scatter: Optional[ScatterDescription]) -> Any:
#     if isinstance(node, StepTagInput):
#         return unwrap(tinput_name, node.source_map[0].source, scatter)
    
#     elif isinstance(node, InputNodeSelector):
#         if channels.exists(janis_uuid=node.input_node.uuid):
#             return unwrap_channel(
#                 tinput_name=tinput_name, 
#                 upstream_node=node.input_node, 
#                 scatter=scatter
#             )
#         elif params.exists(janis_uuid=node.input_node.uuid):
#             param = params.get(janis_uuid=node.input_node.uuid)
#             return [f'params.{param.name}']
    
#     elif isinstance(node, StepOutputSelector):
#         return unwrap_connection(
#             tinput_name=tinput_name, 
#             upstream_step=node.node, 
#             upstream_out=node.tag,
#             scatter=scatter
#         )
    
#     elif isinstance(node, IndexOperator):
#         args = unwrap(tinput_name, node.args[0], scatter)
#         index = node.args[1]
#         args = [f"{arg}[{index}]" for arg in args]
#         return args
    
#     elif isinstance(node, FirstOperator):
#         # TODO implement properly
#         return unwrap(tinput_name, node.args[0][0], scatter)
    
#     elif isinstance(node, AliasSelector):
#         return unwrap(tinput_name, node.inner_selector, scatter)
    
#     elif isinstance(node, list):
#         # TODO this is first selector?
#         if len(node) > 0:
#             return node[0]
#         return None
#     else:
#         raise NotImplementedError




# def unwrap_channel(tinput_name: str, upstream_node: InputNode, scatter: Optional[ScatterDescription]) -> list[str]:
#     """
#     ch_name                     = same type (most cases)
#     ch_name.collect()           = singles -> array (workflow input array channel creation)
#     ch_name.flatten()           = array -> singles (scatter.dot)
#     cartesian_cross.ch_subname  = scatter.cross  
#     """
#     args: list[str] = []
#     relevant_channels = channels.getall(janis_uuid=upstream_node.uuid)
    
#     # arrays of secondaries
#     if relevant_channels and len(relevant_channels) > 1:
#         for ch in relevant_channels:
#             arg = get_channel_expression(
#                 tinput_name=tinput_name, 
#                 channel_name=ch.name,
#                 upstream_dtype=upstream_node.datatype,
#                 scatter=scatter
#             )
#             args.append(arg)
    
#     # everything else
#     elif relevant_channels and len(relevant_channels) == 1:
#         arg = get_channel_expression(
#             tinput_name=tinput_name, 
#             channel_name=relevant_channels[0].name,
#             upstream_dtype=upstream_node.datatype,
#             scatter=scatter
#         )
#         args.append(arg)
#     return args


# def unwrap_connection(tinput_name: str, upstream_step: StepNode, upstream_out: str, scatter: Optional[ScatterDescription]) -> list[str]:
#     # if scatter & output type is Array, use .flatten()
#     args: list[str] = []
#     conn_out = [x for x in upstream_step.tool.tool_outputs() if x.tag == upstream_out][0]

#     # arrays of secondaries
#     if nfgen_utils.is_array_secondary_type(conn_out.outtype):
#         raise NotImplementedError
#     else:
#         upstream_step_id = to_case(upstream_step.id(), settings.NEXTFLOW_PROCESS_CASE)
#         channel_name: str = f'{upstream_step_id}.out.{upstream_out}'
#         arg = get_channel_expression(
#                 tinput_name=tinput_name, 
#                 channel_name=channel_name,
#                 upstream_dtype=conn_out.outtype,
#                 scatter=scatter
#             )
#         args.append(arg)
#     return args

# def get_channel_expression(
#     tinput_name: str, 
#     channel_name: str,
#     upstream_dtype: DataType, 
#     scatter: Optional[ScatterDescription]
#     ) -> str:
#     # scatter
#     if scatter and tinput_name in scatter.fields:
#         # ch_bams -> ch_cartesian_cross.bams
#         if scatter.method == ScatterMethod.cross:
#             return cartesian_cross_subname(channel_name)
#         # ch_bams -> ch_bams.flatten()
#         elif scatter.method == ScatterMethod.dot:
#             if upstream_dtype.is_array():
#                 return f'{channel_name}.flatten()'
#     # everything else
#     return channel_name




"""

# def unwrap_source(src: StepTagInput, scatter: Optional[ScatterDescription]) -> list[str]:
#     source = src.source_map[0].source
#     if isinstance(source, InputNodeSelector):
#         if channels.exists(janis_uuid=source.input_node.uuid):
#             return unwrap_channel(
#                 tinput_name=src.ftag, 
#                 upstream_node=source.input_node, 
#                 scatter=scatter
#             )
#         elif params.exists(janis_uuid=source.input_node.uuid):
#             param = params.get(janis_uuid=source.input_node.uuid)
#             return [f'params.{param.name}']
    
#     elif isinstance(source, StepOutputSelector):
#         return unwrap_connection(
#             tinput_name=src.ftag, 
#             upstream_step=source.node, 
#             upstream_out=source.tag,
#             scatter=scatter
#         )
    
#     elif isinstance(source, AliasSelector):
#         return unwrap_connection(
#             tinput_name=src.ftag,
#             upstream_step=source.inner_selector.node, 
#             upstream_out=source.inner_selector.tag,
#             scatter=scatter
#         )
#     elif isinstance(source, FirstOperator):
#         return unwrap_connection(
#             tinput_name=src.ftag,
#             upstream_step=source.inner_selector.node, 
#             upstream_out=source.inner_selector.tag,
#             scatter=scatter
#         )
    
#     # IndexOperators are annoying
#     elif isinstance(source, IndexOperator):
#         if isinstance(source.args[0], StepOutputSelector):
#             args = unwrap_connection(
#                 tinput_name=src.ftag, 
#                 upstream_step=source.args[0].node, 
#                 upstream_out=source.args[0].tag,
#                 scatter=scatter
#             )
#         elif isinstance(source.args[0], InputNodeSelector):
#             args = unwrap_channel(
#                 tinput_name=src.ftag, 
#                 upstream_node=source.args[0].input_node, 
#                 scatter=scatter
#             )
#         else:
#             raise NotImplementedError
#         index = source.args[1]
#         args = [f"{arg}[{index}]" for arg in args]
#         return args

#     else:
#         raise NotImplementedError


# def unwrap_param() -> list[str]:

"""