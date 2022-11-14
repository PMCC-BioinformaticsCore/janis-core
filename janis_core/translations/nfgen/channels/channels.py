

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import Workflow
from janis_core.types import Array, File

from ..common import NFBase
from ..casefmt import to_case
from .. import utils
from .. import settings
from .. import params
from ..params import Param


### factory 

def register(workflow: Workflow) -> None:
    wfinp_ids = utils.get_channel_input_ids(workflow)
    wfinps = utils.items_with_id(list(workflow.input_nodes.values()), wfinp_ids)
    for inp in wfinps:
        # array of secondaries
        # split to channel per associated param
        # for example - Array(indexedBam) will have 2 params, and 2 channels:
        # params.indexed_bam_bams -> ch_indexed_bam_bams
        # params.indexed_bam_bais -> ch_indexed_bam_bais
        if isinstance(inp.datatype, Array):
            subtype = inp.datatype.subtype()
            if isinstance(subtype, File) and subtype.has_secondary_files():
                register_channels_secondaries_array(inp)
            else:
                register_channel(inp)
        
        # other File type workflow inputs get single channel
        else:
            register_channel(inp)
    
def register_channels_secondaries_array(inp: InputNode) -> None:
    subtype = inp.datatype.subtype()
    exts: list[str] = []
    exts.append(subtype.extension)
    exts += subtype.secondary_files()
    exts = [x.split('.')[-1] for x in exts]
    for ext in exts:
        param_name = f'{inp.id()}_{ext}s'
        param_name = to_case(param_name, settings.NEXTFLOW_PARAM_CASE)
        param = params.get(param_name)
        add(
            varname=param_name,
            params=[param],
            method=get_channel_method(inp),
            collect=should_collect(inp),
            allow_null=should_allow_null(inp),
            reference=inp.id()
        )
        
def register_channel(inp: InputNode) -> None:
    add(
        varname=inp.id(),
        params=params.getall(reference=inp.id()),
        method=get_channel_method(inp),
        collect=should_collect(inp),
        allow_null=should_allow_null(inp),
        reference=inp.id()
    )

def get_channel_method(wfinp: InputNode) -> str:
    # if utils.is_path(wfinp) and utils.is_file_pair(wfinp):
    #     method = 'fromFilePairs'
    if utils.is_path(wfinp):
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

### main class 

@dataclass
class Channel(NFBase):
    varname: str
    params: list[Param]
    method: str
    collect: bool=False
    allow_null: bool=False
    reference: Optional[str]=None
    condensed: bool=True  
    
    @property
    def name(self) -> str:
        name = to_case(self.varname, settings.NEXTFLOW_CHANNEL_CASE)
        return f'ch_{name}'

    @property
    def source(self) -> str:
        param_names = [f'params.{p.name}' for p in self.params]
        if len(param_names) == 1:
            return param_names[0]
        else:
            return repr(param_names)
    
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

def add(
    varname: str,
    params: list[Param],
    method: str,
    collect: bool,
    allow_null: bool,
    reference: Optional[str]=None,
    ) -> None:
    global channel_register
    new_ch = Channel(varname, params, method, collect, allow_null, reference)
    channel_register.channels[new_ch.name] = new_ch

def exists(name: str) -> bool:
    global channel_register 
    if name in channel_register.channels:
        return True
    return False

def get(name: str) -> Channel:
    global channel_register 
    return channel_register.channels[name]
        
def getall() -> list[Channel]:
    global channel_register 
    return channel_register.ordered_channels

def getstr() -> str:
    global channel_register 
    return channel_register.get_string()

def clear() -> None:
    global channel_register 
    channel_register = ChannelRegister()

