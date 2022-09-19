

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import ToolInput
from janis_core.workflow.workflow import InputNode

from .common import NFBase
from . import utils



def channel_factory(task_input: ToolInput | InputNode) -> ChannelDeclaration:
    name = task_input.id()
    method = get_channel_method(task_input)
    source = get_channel_source(task_input, method)
    collect = should_collect(task_input)
    return ChannelDeclaration(name, method, source, collect)

def get_channel_method(task_input: ToolInput | InputNode) -> str:
    # consumed once or multiple times? (only valid for InputNode)
    if isinstance(task_input, InputNode) and len(utils.get_references(task_input)) > 1:
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


@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        ...

@dataclass
class QueueChannelPriority(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        channels.sort(key=lambda x: x.method != 'value', reverse=True) 
        return channels

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        channels.sort(key=lambda x: x.method == 'fromPath' or x.method == 'fromFilePairs', reverse=True) 
        return channels
    
@dataclass
class Alphabetical(OrderingMethod):
    def order(self, channels: list[ChannelDeclaration]) -> list[ChannelDeclaration]:
        channels.sort(key=lambda x: x.name) 
        return channels

orderers: list[OrderingMethod] = [
    Alphabetical(),
    QueueChannelPriority(),
    FileTypePriority(),
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
        width_col_1 = max([c.width for c in self.channels])
        outstr = ''
        for c in self.ordered_channels:
            outstr += f'{c.name:<{width_col_1}} = {c.get_string()}\n'
        return outstr



class ChannelDeclaration(NFBase):
    condensed: bool = True   # class method - shared by each instance. 

    def __init__(self, name: str, method: str, source: str, collect: bool=False):
        self._name = name
        self.method = method
        self.source = source
        self.collect = collect
    
    @property
    def name(self) -> str:
        return f'ch_{self._name}'
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        return self.get_string_condensed()
        
    def get_string_condensed(self) -> str:
        collect = '.collect()' if self.collect else ''
        return f'Channel.{self.method}( {self.source} ){collect}'

    def get_string_expanded(self) -> str:
        channel_str = ''
        channel_str += 'Channel\n'
        channel_str += f'  .{self.method}( {self.source} )\n'
        channel_str += f'  .collect()\n' if self.collect else ''
        channel_str += f'  .set{{ {self.name} }}\n'
        return channel_str
    





