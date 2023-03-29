

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional
from collections import defaultdict
from copy import deepcopy

from janis_core.types import DataType
from uuid import uuid4


from ..scope import Scope
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
        self.data_structure: dict[str, list[Channel]] = defaultdict(list)

    @property
    def all_channels(self) -> list[Channel]:
        channels: list[Channel] = []
        for scoped_channels in self.data_structure.values():
            channels += scoped_channels
        return order(channels)
    
    def scoped_channels(self, scope: Scope, for_parent: bool=False) -> list[Channel]:
        if for_parent:
            scope_copy = deepcopy(scope)
            scope_copy.items = scope_copy.items[:-1]
            label = scope_copy.to_string()
        else:
            label = scope.to_string()
        channels = self.data_structure[label]
        return order(channels)
    
    def update(self, scope: Scope, channel: Channel) -> None:
        scope_label = scope.to_string()
        self.data_structure[scope_label].append(channel)

    def get_string(self, scope: Scope) -> str:
        is_sub_wf = True if len(scope.labels) > 1 else False
        label = scope.to_string()
        scoped_channels = self.data_structure[label]
        if is_sub_wf:
            return self._render_subwf(scoped_channels)
        else:
            return self._render_main(scoped_channels)

    def _render_main(self, channels: list[Channel]) -> str:
        outstr = ''
        channels = order(channels)
        channels = [ch for ch in channels if ch.define]
        width_col_1 = max([ch.width for ch in channels])
        for ch in channels:
            outstr += f'{ch.name:<{width_col_1}} = {ch.get_string()}\n'
        return outstr

    def _render_subwf(self, channels: list[Channel]) -> str:
        outstr = ''
        channels = order(channels)
        for ch in channels:
            outstr += f'{ch.name}\n'
        return outstr


channel_register = ChannelRegister()

def add(
    scope: Scope,
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
    channel_register.update(scope, new_ch)

def exists(scope: Scope, janis_uuid: str, for_parent: bool=False) -> bool:
    global channel_register
    scoped_channels = channel_register.scoped_channels(scope, for_parent=for_parent)
    for ch in scoped_channels:
        if ch.janis_uuid == janis_uuid:
            return True
    return False

def get(scope: Scope, janis_uuid: str, for_parent: bool=False) -> Channel:
    global channel_register
    scoped_channels = channel_register.scoped_channels(scope, for_parent=for_parent)
    for ch in scoped_channels:
        if ch.janis_uuid == janis_uuid:
            return ch
    # exception - channel not found. 
    label = scope.to_string()
    raise LookupError(f'no channel for janis_uuid within scope {label}. use .exists() to check')

def getall(scope: Scope) -> list[Channel]:
    global channel_register
    return channel_register.scoped_channels(scope)

def getstr(scope: Scope) -> str:
    global channel_register 
    return channel_register.get_string(scope)

def clear() -> None:
    global channel_register 
    channel_register = ChannelRegister()





