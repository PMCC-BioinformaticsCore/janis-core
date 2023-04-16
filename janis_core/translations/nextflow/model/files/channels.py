

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from abc import ABC, abstractmethod
from dataclasses import dataclass


### main class 

@dataclass
class NFChannelDefinition:
    name: str
    source: str
    method: str
    operations: Optional[str]=None

    @property
    def name_width(self) -> int:
        return len(self.name)

    def get_string(self) -> str:
        return f'Channel.{self.method}( {self.source} ){self.operations}'
    

@dataclass
class NFChannelDefinitionBlock:
    channels: list[NFChannelDefinition]

    @property
    def def_width(self) -> int:
        return max(ch.name_width for ch in self.channels) + 1

    def get_string(self) -> str:
        out: str = ''
        ordered_channels = order_nf_channels(self.channels)
        for ch in ordered_channels:
            line = f'{ch.name:<{self.def_width}} = {ch.get_string()}\n'
            out += line
        return out




### ordering 

@dataclass
class OrderingMethod(ABC):
    @abstractmethod
    def order(self, channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
        ...

@dataclass
class QueueChannelPriority(OrderingMethod):
    def order(self, channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
        return sorted(channels, key=lambda x: x.method != 'value', reverse=True) 

@dataclass
class FileTypePriority(OrderingMethod):
    def order(self, channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
        return sorted(channels, key=lambda x: x.method == 'fromPath' or x.method == 'fromFilePairs', reverse=True) 

@dataclass
class MandatoryPriority(OrderingMethod):
    def order(self, channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
        return sorted(channels, key=lambda x: 'ifEmpty' in x.operations if x.operations else False, reverse=True)
    
@dataclass
class Alphabetical(OrderingMethod):
    def order(self, channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
        return sorted(channels, key=lambda x: x.name) 

orderers: list[OrderingMethod] = [
    Alphabetical(),
    QueueChannelPriority(),
    FileTypePriority(),
    MandatoryPriority()
]

def order_nf_channels(channels: list[NFChannelDefinition]) -> list[NFChannelDefinition]:
    for orderer in orderers:
        channels = orderer.order(channels)
    return channels


