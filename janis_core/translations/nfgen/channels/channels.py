

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import Workflow

from ..common import NFBase
from .. import utils


### factory 

def register(workflow: Workflow) -> None:
    for wfinp in utils.get_workflow_inputs(workflow):
        add(wfinp)

def get_channel_method(wfinp: InputNode) -> str:
    if utils.is_path(wfinp) and utils.is_file_pair(wfinp):
        method = 'fromFilePairs'
    elif utils.is_path(wfinp):
        method = 'fromPath'
    else: 
        method = 'of'
    return method

def get_channel_source(wfinp: InputNode) -> str:
    return f'params.{wfinp.id()}'

def should_collect(wfinp: InputNode) -> bool:
    if wfinp.datatype.is_array():
        return True
    elif wfinp.default is not None:
        if isinstance(wfinp.default, list):
            return True
    return False

def should_allow_null(wfinp: InputNode) -> bool:
    if wfinp.datatype.optional:
        return True
    return False


### ordering 

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, channels: list[Channel]) -> list[Channel]:
        ...

@dataclass
class QueueChannelPriority(OrderingMethod):
    def order(self, channels: list[Channel]) -> list[Channel]:
        return sorted(channels, key=lambda x: x.method != 'value', reverse=True) 

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, channels: list[Channel]) -> list[Channel]:
        return sorted(channels, key=lambda x: x.method == 'fromPath' or x.method == 'fromFilePairs', reverse=True) 

@dataclass
class MandatoryPriority(OrderingMethod):
    def order(self, channels: list[Channel]) -> list[Channel]:
        return sorted(channels, key=lambda x: x.allow_null == False, reverse=True)
    
@dataclass
class Alphabetical(OrderingMethod):
    def order(self, channels: list[Channel]) -> list[Channel]:
        return sorted(channels, key=lambda x: x.name) 

orderers: list[OrderingMethod] = [
    Alphabetical(),
    QueueChannelPriority(),
    FileTypePriority(),
    MandatoryPriority()
]

### main class 

@dataclass
class Channel(NFBase):
    wfinp_name: str
    method: str
    source: str
    collect: bool=False
    allow_null: bool=False
    condensed: bool=True  
    
    @property
    def name(self) -> str:
        return f'ch_{self.wfinp_name}'
    
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
    

### register for channels

class ChannelRegister(NFBase):
    def __init__(self):
        self.channels: dict[str, Channel] = {}

    @property
    def ordered_channels(self) -> list[Channel]:
        channels = list(self.channels.values())
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


channel_register = ChannelRegister()

def add(wfinp: InputNode) -> None:
    global channel_register 
    new_ch = Channel(
            wfinp_name=wfinp.id(),
            method = get_channel_method(wfinp),
            source = get_channel_source(wfinp),
            collect = should_collect(wfinp),
            allow_null = should_allow_null(wfinp),
        )
    channel_register.channels[new_ch.wfinp_name] = new_ch

def exists(wfinp_name: str) -> bool:
    global channel_register 
    if wfinp_name in channel_register.channels:
        return True
    return False

def get(wfinp_name: str) -> Channel:
    global channel_register 
    return channel_register.channels[wfinp_name]
        
def getall() -> list[Channel]:
    global channel_register 
    return channel_register.ordered_channels

def getstr() -> str:
    global channel_register 
    return channel_register.get_string()


