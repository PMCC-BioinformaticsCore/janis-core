


from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import DataType
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType


@dataclass
class NFProcessInput(ABC):
    name: str           # the 
    tinput_id: str      # the janis tool input id
    dtype: DataType
    
    @abstractmethod
    def get_string(self) -> str:
        ...



@dataclass
class NFValProcessInput(NFProcessInput):

    def get_string(self) -> str:
        return f'val {self.name}'


@dataclass
class NFPathProcessInput(NFProcessInput):
    presents_as: Optional[str]=None

    @property
    def stage_as(self) -> str:
        # presents_as takes precedent? this might be really bad.
        dtt = utils.get_dtt(self.dtype)
        if self.presents_as:
            expr = f", stageAs: '{self.presents_as}'"
        elif dtt == DTypeType.FILE_ARRAY and self.dtype.optional:
            expr = f", stageAs: '{self.name}??/*'"
        elif dtt == DTypeType.SECONDARY and self.dtype.optional:
            expr = f", stageAs: '{self.name}/*'"
        elif dtt == DTypeType.FILE and self.dtype.optional:
            expr = f", stageAs: '{self.name}/*'"
        else:
            expr = ''
        return expr
    
    def get_string(self) -> str:
        return f'path {self.name}{self.stage_as}'


@dataclass
class NFScriptProcessInput(NFProcessInput):
    presents_as: str

    def get_string(self) -> str:
        return f"path {self.name}, stageAs: '{self.presents_as}'"

@dataclass
class NFPythonToolProcessInput(NFProcessInput):

    def get_string(self) -> str:
        return f'path {self.name}'
    

@dataclass
class NFTupleProcessInput(NFProcessInput):
    subnames: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for subname in self.subnames:
            out += f'path({subname}{self.stage_as(subname)}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out
    
    def stage_as(self, subname: str) -> str:
        # optional secondary arrays
        expr = ''
        # if utils.is_secondary_type(self.dtype) and self.dtype.optional:
        #     expr = f", stageAs: '{subname}/*'"

        # # optional file pair arrays
        # elif utils.is_file_pair_type(self.dtype) and self.dtype.optional:
        #     expr = f", stageAs: '{subname}/*'"

        return expr

    def get_string(self) -> str:
        return f'tuple {self.fields}'