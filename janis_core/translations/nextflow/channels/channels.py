

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from janis_core.types import DataType
from uuid import uuid4

from .. import naming


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
        return sorted(channels, key=lambda x: 'ifEmpty' in x.operations if x.operations else False, reverse=True)
    
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

def order(channels: list[Channel]) -> list[Channel]:
    for orderer in orderers:
        channels = orderer.order(channels)
    return channels

### main class 

@dataclass
class Channel:
    name: str
    source: str
    method: str
    operations: Optional[str]=None
    janis_uuid: Optional[str]=None
    define: bool=False

    def __post_init__(self):
        self.uuid = uuid4() 
    
    @property
    def definition(self) -> str:
        return f'{self.name} = {self.get_string()}'

    # @property
    # def source(self) -> str:
    #     param_names = [f'params.{p.name}' for p in self.params]
    #     if len(param_names) == 1:
    #         return param_names[0]
    #     elif len(param_names) == 2:
    #         return f"[{', '.join(param_names)}]"
    #     elif len(param_names) >= 3:
    #         expr = ''
    #         expr += '[\n'
    #         for pname in param_names:
    #             expr += f'{settings.NF_INDENT}{pname},\n'
    #         expr += ']'
    #         return expr
    #     else:
    #         raise RuntimeError('DEV: should have 1+ param but didnt.')
    
    @property
    def width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        return f'Channel.{self.method}( {self.source} ){self.operations}'

    

### register for channels

class ChannelRegister:
    def __init__(self):
        self.channels: list[Channel] = []

    @property
    def ordered_channels(self) -> list[Channel]:
        return order(self.channels)

    def get_string(self) -> str:
        outstr = ''
        channels = self.ordered_channels
        channels = [ch for ch in channels if ch.define]
        width_col_1 = max([ch.width for ch in channels])
        for ch in channels:
            outstr += f'{ch.name:<{width_col_1}} = {ch.get_string()}\n'
        return outstr


channel_register = ChannelRegister()

def add(
    janis_tag: str,
    method: str,
    source: str,
    operations: Optional[str]=None,
    name_override: Optional[str]=None,
    janis_dtype: Optional[DataType]=None,
    janis_uuid: Optional[str]=None,
    define: bool=False
    ) -> None:
    global channel_register
    # channel name
    name = naming.constructs.gen_varname_channel(janis_tag, name_override, janis_dtype)
    # create channel
    new_ch = Channel(name, source, method, operations, janis_uuid, define)
    # add channel
    channel_register.channels.append(new_ch)

def exists(janis_uuid: str) -> bool:
    global channel_register
    channels = getall(janis_uuid)
    if channels:
        return True
    return False

def get(janis_uuid: str) -> Channel:
    global channel_register 
    return getall(janis_uuid)[0]

def getall(janis_uuid: Optional[str]=None) -> list[Channel]:
    global channel_register
    channels = channel_register.ordered_channels
    if janis_uuid:
        channels = [x for x in channels if x.janis_uuid == janis_uuid]
    return channels

def getstr() -> str:
    global channel_register 
    return channel_register.get_string()

def clear() -> None:
    global channel_register 
    channel_register = ChannelRegister()





