


from typing import Optional

from janis_core import TInput, Workflow
from janis_core.types import File, Filename, Directory, DataType

from janis_core import settings
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ...model.files.channels import NFChannelDefinition
from ...model.files.channels import NFChannelDefinitionBlock
from ...model.workflow import NFWorkflow
from ...model.workflow import NFMainWorkflow

from ... import naming 
from ... import task_inputs
from ...task_inputs import TaskInputType

INDENT = settings.translate.nextflow.NF_INDENT
CROSS_CHANNEL_NAME = 'ch_cartesian_cross'



def gen_channels_block(nf_workflow: NFWorkflow, wf: Workflow) -> Optional[NFChannelDefinitionBlock]:
    # for each param
    # if it qualifies to create a channel, create channel declaration 
    
    channel_definitions: list[NFChannelDefinition] = []
    channel_block: Optional[NFChannelDefinitionBlock] = None
    
    if isinstance(nf_workflow, NFMainWorkflow):
        for input_node in wf.tool_inputs():
            if should_create_channel_definition(input_node, wf):
                generator = ChannelDefinitionGenerator(input_node, wf)
                ch_def = generator.register()
                channel_definitions.append(ch_def)

    if channel_definitions:
        channel_block = NFChannelDefinitionBlock(channel_definitions)
    
    return channel_block

def should_create_channel_definition(input_node: TInput, wf: Workflow) -> bool:
    if not task_inputs.exists(wf.id(), input_node):
        return False
    
    task_input = task_inputs.get(wf.id(), input_node)
    if task_input.ti_type in (TaskInputType.STATIC, TaskInputType.IGNORED, TaskInputType.LOCAL):
        return False
    
    dtt = utils.get_dtt(input_node.intype)
    if dtt not in [
        DTypeType.SECONDARY_ARRAY,
        DTypeType.SECONDARY,
        DTypeType.FILE_PAIR_ARRAY,
        DTypeType.FILE_PAIR,
        DTypeType.FILE_ARRAY,
        DTypeType.FILE,
    ]:
        return False

    if input_node.intype.optional:
        return False
    
    return True


class ChannelDefinitionGenerator:
    def __init__(self, tinput: TInput, wf: Workflow) -> None:
        self.tinput = tinput
        self.wf = wf

    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.tinput.intype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_channel(self.tinput.id())
    
    @property
    def method(self) -> str:
        # if isinstance(self.basetype, (File, Filename, Directory)) and not self.tinput.intype.optional:
        if isinstance(self.basetype, (File, Filename, Directory)):
            return 'fromPath'
        return 'of'

    @property
    def source(self) -> str:
        task_input = task_inputs.get(self.wf.id(), self.tinput)
        param_name = task_input.value
        # @secondaryarrays
        if utils.is_secondary_array_type(self.tinput.intype):
            src = f'{param_name}.flatten()'
        else:
            src = f'{param_name}'
        return src

    def register(self) -> NFChannelDefinition:
        operations = self.get_operations()
        return NFChannelDefinition(
            name=self.name,
            source=self.source,
            method=self.method,
            operations=operations
        )
    
    def get_operations(self) -> str:
        # secondary array
        if utils.is_secondary_array_type(self.tinput.intype):
            ops = self.get_operations_secondary_array()
        
        # secondary
        elif utils.is_secondary_type(self.tinput.intype):
            ops = self.get_operations_secondary()
        
        # file array
        elif isinstance(self.basetype, (File, Filename, Directory)) and self.tinput.intype.is_array():
            ops = self.get_operations_file_array()
        
        # nonfile array
        elif self.tinput.intype.is_array():
            ops = self.get_operations_nonfile_array()
        
        # anything else
        else:
            ops = self.get_operations_generic()
        return ops

    def get_operations_secondary_array(self) -> str:
        exts = utils.get_extensions(self.basetype)
        size = len(exts)
        
        ops: str = ''
        ops += f'.collate( {size} )'
        # if self.tinput.intype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_secondary(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.tinput.intype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_file_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.tinput.intype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_nonfile_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.tinput.intype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_generic(self) -> str:
        ops: str = ''
        # if self.tinput.intype.optional:
        #     ops += '.ifEmpty( null )'
        return ops
    



### DEPRECATED ###

# class ChannelOperation(ABC):
#     pass


# @dataclass 
# class CartesianCrossOperation(ChannelOperation):
#     channel_names: list[str]
#     channel_dtypes: list[DataType]

#     def get_string(self) -> str:
#         lines: list[str] = []

#         # first channel to combine
#         if self.channel_dtypes[0].is_array():
#             lines.append(f'{self.channel_names[0]}.flatten()')  # ch_bams.flatten()
#         else:
#             lines.append(self.channel_names[0])                 # ch_bams
        
#         # other channels to combine
#         for ch_name, ch_dtype in zip(self.channel_names[1:], self.channel_dtypes[1:]):
#             if ch_dtype.is_array():
#                 lines.append(f'.combine({ch_name}).flatten()')  # .combine(ch_bais).flatten()
#             else:
#                 lines.append(f'.combine({ch_name})')            # .combine(ch_bais)
        
#         # separation back to individual channels
#         lines.append('.multiMap { it ->')                       # .multiMap { it ->
#         for i, ch_name in enumerate(self.channel_names):
#             ch_subname = ch_name.replace('ch_', '')             # ch_bams -> bams
#             ch_subname = ch_subname.split('.')[-1]              # STP1.out.bams -> bams
#             lines.append(f'{INDENT}{ch_subname}: it[{i}]')      # bams: it[0] \n bais: it[1]
#         lines.append('}')

#         # set new channel names
#         lines.append(f'.set {{ {CROSS_CHANNEL_NAME} }}')
#         text = '\n'.join(lines)
#         return text

# def gen_scatter_cross_operation(sources: dict[str, Any], scatter: ScatterDescription) -> ChannelOperation:
#     channel_names: list[str] = []
#     channel_dtypes: list[DataType] = []
#     for tinput_name, src in sources.items():
#         if tinput_name in scatter.fields:
#             channel_names.append(resolve_channel_name(src))
#             channel_dtypes.append(resolve_channel_dtype(src))

#     return CartesianCrossOperation(
#         channel_names=channel_names,
#         channel_dtypes=channel_dtypes,
#     )

# def resolve_channel_name(src: StepTagInput) -> str:
#     # workflow input
#     source = src.source_map[0].source
#     if isinstance(source, InputNodeSelector):
#         relevant_channels = getall(janis_uuid=source.input_node.uuid)
#         assert(relevant_channels)
#         assert(len(relevant_channels) == 1)
#         return relevant_channels[0].name
        
#     # step output
#     elif isinstance(source, StepOutputSelector):
#         conn_step: StepNode     = source.node
#         conn_step_id: str       = naming.constructs.gen_varname_process(conn_step.id())
#         conn_out_id: str        = source.tag
#         channel_name: str       = f'{conn_step_id}.out.{conn_out_id}'
#         return channel_name
    
#     raise NotImplementedError

# def resolve_channel_dtype(src: StepTagInput) -> DataType:
#     # workflow input
#     source = src.source_map[0].source
#     if isinstance(source, InputNodeSelector):
#         return source.input_node.datatype
#     # step output
#     elif isinstance(source, StepOutputSelector):
#         conn_step: StepNode     = source.node
#         conn_out_id: str        = source.tag
#         conn_out: ToolOutput    = [x for x in conn_step.tool.outputs() if x.tag == conn_out_id][0]
#         return conn_out.output_type
    
#     raise NotImplementedError