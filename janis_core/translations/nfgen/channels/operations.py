


from typing import Any
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import ToolOutput
from janis_core.workflow.workflow import StepNode
from janis_core.types import DataType
from janis_core.utils.scatter import ScatterDescription
from janis_core.graph.steptaginput import StepTagInput
from janis_core.operators.selectors import InputNodeSelector
from janis_core.operators.selectors import StepOutputSelector
from ..casefmt import to_case
from ..common import NFBase
from .. import settings
from .channels import getall

INDENT = settings.NEXTFLOW_INDENT
CROSS_CHANNEL_NAME = 'ch_cartesian_cross'




class ChannelOperation(NFBase):
    
    @abstractmethod
    def get_string(self) -> str:
        ...

@dataclass 
class CartesianCrossOperation(ChannelOperation):
    channel_names: list[str]
    channel_dtypes: list[DataType]

    def get_string(self) -> str:
        lines: list[str] = []

        # first channel to combine
        if self.channel_dtypes[0].is_array():
            lines.append(f'{self.channel_names[0]}.flatten()')  # ch_bams.flatten()
        else:
            lines.append(self.channel_names[0])                 # ch_bams
        
        # other channels to combine
        for ch_name, ch_dtype in zip(self.channel_names[1:], self.channel_dtypes[1:]):
            if ch_dtype.is_array():
                lines.append(f'.combine({ch_name}).flatten()')  # .combine(ch_bais).flatten()
            else:
                lines.append(f'.combine({ch_name})')            # .combine(ch_bais)
        
        # separation back to individual channels
        lines.append('.multiMap { it ->')                       # .multiMap { it ->
        for i, ch_name in enumerate(self.channel_names):
            ch_subname = ch_name.replace('ch_', '')             # ch_bams -> bams
            ch_subname = ch_subname.split('.')[-1]              # STP1.out.bams -> bams
            lines.append(f'{INDENT}{ch_subname}: it[{i}]')      # bams: it[0] \n bais: it[1]
        lines.append('}')

        # set new channel names
        lines.append(f'.set {{ {CROSS_CHANNEL_NAME} }}')
        text = '\n'.join(lines)
        return text




def gen_scatter_cross_operation(sources: dict[str, Any], scatter: ScatterDescription) -> ChannelOperation:
    channel_names: list[str] = []
    channel_dtypes: list[DataType] = []
    for tinput_name, src in sources.items():
        if tinput_name in scatter.fields:
            channel_names.append(resolve_channel_name(src))
            channel_dtypes.append(resolve_channel_dtype(src))

    return CartesianCrossOperation(
        channel_names=channel_names,
        channel_dtypes=channel_dtypes,
    )

def resolve_channel_name(src: StepTagInput) -> str:
    # workflow input
    source = src.source_map[0].source
    if isinstance(source, InputNodeSelector):
        relevant_channels = getall(janis_uuid=source.input_node.uuid)
        if relevant_channels and len(relevant_channels) == 1:
            return relevant_channels[0].name
        else:
            raise NotImplementedError
    # step output
    elif isinstance(source, StepOutputSelector):
        conn_step: StepNode     = source.node
        conn_step_id: str       = to_case(conn_step.id(), settings.NEXTFLOW_PROCESS_CASE)
        conn_out_id: str        = source.tag
        channel_name: str       = f'{conn_step_id}.out.{conn_out_id}'
        return channel_name
    else:
        raise NotImplementedError

def resolve_channel_dtype(src: StepTagInput) -> DataType:
    # workflow input
    source = src.source_map[0].source
    if isinstance(source, InputNodeSelector):
        return source.input_node.datatype
    # step output
    elif isinstance(source, StepOutputSelector):
        conn_step: StepNode     = source.node
        conn_out_id: str        = source.tag
        conn_out: ToolOutput    = [x for x in conn_step.tool.outputs() if x.tag == conn_out_id][0]
        return conn_out.output_type
    else:
        raise NotImplementedError