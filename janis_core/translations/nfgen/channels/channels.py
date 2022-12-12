

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from janis_core.workflow.workflow import InputNode
from janis_core.types import File
from uuid import uuid4

from ..common import NFBase
from ..casefmt import to_case
from ..params import Param
from .. import nfgen_utils
from .. import settings


def get_channel_method(wfinp: InputNode) -> str:
    # if utils.is_path(wfinp) and utils.is_file_pair(wfinp):
    #     method = 'fromFilePairs'
    if nfgen_utils.is_path(wfinp):
        method = 'fromPath'
    else: 
        method = 'of'
    return method

def should_collect(wfinp: InputNode) -> bool:
    if wfinp.datatype.is_array():
        return True
    elif wfinp.default is not None:
        if isinstance(wfinp.default, list):
            return True
    elif isinstance(wfinp.datatype, File) and wfinp.datatype.has_secondary_files():
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

def order(channels: list[Channel]) -> list[Channel]:
    for orderer in orderers:
        channels = orderer.order(channels)
    return channels

### main class 

@dataclass
class Channel(NFBase):
    ref_name: str
    ref_scope: list[str]
    params: list[Param]
    method: str
    collect: bool=False
    allow_null: bool=False
    name_override: Optional[str]=None
    janis_uuid: Optional[str]=None
    define: bool=False

    def __post_init__(self):
        self.uuid = uuid4() 
    
    @property
    def name(self) -> str:
        if self.name_override:
            base = self.name_override
        else:
            base = self.ref_name
        base = to_case(base, settings.NF_CHANNEL_CASE)
        full = f'ch_{base}'
        return full

    @property
    def declaration(self) -> str:
        return f'{self.name} = {self.get_string()}'

    @property
    def source(self) -> str:
        param_names = [f'params.{p.name}' for p in self.params]
        if len(param_names) == 1:
            return param_names[0]
        else:
            return ', '.join(param_names)
    
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
    var_name: str,
    var_scope: list[str],
    params: list[Param],
    method: str,
    collect: bool,
    allow_null: bool,
    name_override: Optional[str]=None,
    janis_uuid: Optional[str]=None,
    define: bool=False
    ) -> None:
    global channel_register
    new_ch = Channel(var_name, var_scope, params, method, collect, allow_null, name_override, janis_uuid, define)
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





