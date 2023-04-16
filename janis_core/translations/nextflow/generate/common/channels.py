


# from typing import Any
# from abc import ABC
# from dataclasses import dataclass

# from janis_core import ToolOutput, TInput, Tool
# from janis_core.types import File, Filename, Directory, DataType
# from janis_core.workflow.workflow import InputNode, StepNode
# from janis_core.utils.scatter import ScatterDescription
# from janis_core.graph.steptaginput import StepTagInput
# from janis_core.operators.selectors import InputNodeSelector
# from janis_core.operators.selectors import StepOutputSelector

# from janis_core import settings
# from janis_core import translation_utils as utils

# from ... import naming 
# from ...model.files.channels import NFChannelDefinition
# from ...model.files.channels import NFChannelDefinitionBlock

# INDENT = settings.translate.nextflow.NF_INDENT
# CROSS_CHANNEL_NAME = 'ch_cartesian_cross'




# # def add(
# #     scope: Scope,
# #     janis_tag: str,
# #     method: str,
# #     source: str,
# #     operations: Optional[str]=None,
# #     name_override: Optional[str]=None,
# #     janis_dtype: Optional[DataType]=None,
# #     janis_uuid: Optional[str]=None,
# #     define: bool=False
# #     ) -> None:
# #     global channel_register
# #     # channel name
# #     name = naming.constructs.gen_varname_channel(janis_tag, name_override, janis_dtype)
# #     # create channel
# #     new_ch = ChannelDefinition(name, source, method, operations, janis_uuid, define)
# #     # add channel
# #     channel_register.update(scope, new_ch)



# def create_channel_definitions_block(tool: Tool, tinput_ids: set[str]) -> NFChannelDefinitionBlock:
#     channels: list[NFChannelDefinition] = [] 
#     for tinput_id in tinput_ids:
#         tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
#         ch_def = create_channel_definition(tinput)
#         channels.append(ch_def)
#     return NFChannelDefinitionBlock(channels)

# def create_channel_definition(tinput) -> NFChannelDefinition:
#     pass
    

# class ChannelRegistrationManager:
#     def __init__(self, tinput: TInput) -> None:
#         pass


#     def __init__(self, inp: InputNode, scope: Scope) -> None:
#         self.inp = inp
#         self.scope = scope

#     @property
#     def basetype(self) -> DataType:
#         basetype = utils.get_base_type(self.inp.datatype)
#         basetype = utils.ensure_single_type(basetype)
#         return basetype

#     @property
#     def is_subworkflow(self) -> bool:
#         if self.scope.labels != [settings.translate.nextflow.NF_MAIN_NAME]:
#             return True
#         return False
    
#     @property
#     def method(self) -> str:
#         # if isinstance(self.basetype, (File, Filename, Directory)) and not self.inp.datatype.optional:
#         if isinstance(self.basetype, (File, Filename, Directory)):
#             return 'fromPath'
#         return 'of'

#     @property
#     def source(self) -> str:
#         if self.is_subworkflow:
#             src = ''
#         else:
#             if not params.exists(self.inp.uuid):
#                 print()
#             param_name = params.get(self.inp.uuid).name
#             # @secondaryarrays
#             if utils.is_array_secondary_type(self.inp.datatype):
#                 src = f'params.{param_name}.flatten()'
#             else:
#                 src = f'params.{param_name}'
#         return src

#     def register(self) -> None:
#         operations = self.get_operations()
#         channels.add(
#             scope=self.scope,
#             janis_tag=self.inp.id(),
#             method=self.method,
#             source=self.source,
#             operations=operations,
#             janis_dtype=self.inp.datatype,
#             janis_uuid=self.inp.uuid,
#             define=False if self.is_subworkflow else True
#         )

#     def get_operations(self) -> str:
#         # secondary array
#         if utils.is_array_secondary_type(self.inp.datatype):
#             ops = self.get_operations_secondary_array()
        
#         # secondary
#         elif utils.is_secondary_type(self.inp.datatype):
#             ops = self.get_operations_secondary()
        
#         # file array
#         elif isinstance(self.basetype, (File, Filename, Directory)) and self.inp.datatype.is_array():
#             ops = self.get_operations_file_array()
        
#         # nonfile array
#         elif self.inp.datatype.is_array():
#             ops = self.get_operations_nonfile_array()
        
#         # anything else
#         else:
#             ops = self.get_operations_generic()
#         return ops

#     def get_operations_secondary_array(self) -> str:
#         exts = utils.get_extensions(self.basetype)
#         size = len(exts)
        
#         ops: str = ''
#         ops += f'.collate( {size} )'
#         # if self.inp.datatype.optional:
#         #     ops += '.ifEmpty( null )'
#         return ops

#     def get_operations_secondary(self) -> str:
#         ops: str = ''
#         ops += '.toList()'
#         # if self.inp.datatype.optional:
#         #     ops += '.ifEmpty( null )'
#         return ops

#     def get_operations_file_array(self) -> str:
#         ops: str = ''
#         ops += '.toList()'
#         # if self.inp.datatype.optional:
#         #     ops += '.ifEmpty( null )'
#         return ops

#     def get_operations_nonfile_array(self) -> str:
#         ops: str = ''
#         ops += '.toList()'
#         # if self.inp.datatype.optional:
#         #     ops += '.ifEmpty( null )'
#         return ops

#     def get_operations_generic(self) -> str:
#         ops: str = ''
#         # if self.inp.datatype.optional:
#         #     ops += '.ifEmpty( null )'
#         return ops
    




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