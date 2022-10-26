

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import ToolInput
from janis_core.workflow.workflow import InputNode

from .common import NFBase
from . import utils



def channel_factory(task_input: ToolInput | InputNode) -> ChannelDeclaration:
    wf_inp_name = task_input.id()
    method = get_channel_method(task_input)
    source = get_channel_source(task_input, method)
    collect = should_collect(task_input)
    allow_null = should_allow_null(task_input)
    return ChannelDeclaration(wf_inp_name, method, source, collect, allow_null)

def get_channel_method(task_input: ToolInput | InputNode) -> str:
    # consumed once or multiple times? (only valid for InputNode)
    if isinstance(task_input, InputNode) and len(utils.get_input_references(task_input)) > 1:
        method = 'value'
    elif utils.is_path(task_input) and utils.is_file_pair(task_input):
        method = 'fromFilePairs'
    elif utils.is_path(task_input):
        method = 'fromPath'
    else: 
        method = 'of'
    return method

def get_channel_source(task_input: ToolInput | InputNode, method: str) -> str:
    if method == 'value' and utils.is_path(task_input):
        source = f"file('params.{task_input.id()}')"
    else:
        source = f'params.{task_input.id()}'
    return source

def should_collect(task_input: ToolInput | InputNode) -> bool:
    if task_input.default is not None:
        if isinstance(task_input.default, list):
            return True
    return False

def should_allow_null(task_input: ToolInput | InputNode) -> bool:
    if isinstance(task_input, ToolInput):
        dtype = task_input.input_type
    else:
        dtype = task_input.datatype
    if dtype.optional:
        return True
    return False


@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        ...

@dataclass
class QueueChannelPriority(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        return sorted(channels, key=lambda x: x.method != 'value', reverse=True) 

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        return sorted(channels, key=lambda x: x.method == 'fromPath' or x.method == 'fromFilePairs', reverse=True) 

@dataclass
class MandatoryPriority(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        return sorted(channels, key=lambda x: x.allow_null == False, reverse=True)
    
@dataclass
class Alphabetical(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        return sorted(channels, key=lambda x: x.name) 

orderers: list[OrderingMethod] = [
    Alphabetical(),
    QueueChannelPriority(),
    FileTypePriority(),
    MandatoryPriority()
]


class ChannelDeclarationBlock(NFBase):
    def __init__(self, channels: list[ChannelDeclaration]):
        self.channels = channels

    @property
    def ordered_channels(self) -> list[ChannelDeclaration]:
        channels = self.channels
        for orderer in orderers:
            channels = orderer.order(channels)
        return channels

    def get_string(self) -> str:
        outstr = ''
        channels = self.ordered_channels
        width_col_1 = max([c.width for c in channels])
        for c in channels:
            outstr += f'{c.name:<{width_col_1}} = {c.get_string()}\n'
        return outstr



class ChannelDeclaration(NFBase):
    condensed: bool = True   # class method - shared by each instance. 

    def __init__(self, wf_inp_name: str, method: str, source: str, collect: bool=False, allow_null: bool=False):
        self.wf_inp_name = wf_inp_name
        self.method = method
        self.source = source
        self.collect = collect
        self.allow_null = allow_null
    
    @property
    def name(self) -> str:
        return f'ch_{self.wf_inp_name}'
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        return self.get_string_condensed()
        
    def get_string_condensed(self) -> str:
        collect = '.collect()' if self.collect else ''
        ifempty = '.ifEmpty( null )' if self.allow_null else ''
        return f'Channel.{self.method}( {self.source} ){collect}{ifempty}'

    def get_string_expanded(self) -> str:
        channel_str = ''
        channel_str += 'Channel\n'
        channel_str += f'  .{self.method}( {self.source} )\n'
        channel_str += f'  .collect()\n' if self.collect else ''
        channel_str += f'  .ifEmpty( null )\n' if self.allow_null else ''
        channel_str += f'  .set{{ {self.name} }}\n'
        return channel_str
    





