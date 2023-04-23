


from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import DataType, File, Directory
from janis_core import translation_utils as utils


@dataclass
class NFProcessInput(ABC):
    name: str
    tinput_id: str
    dtype: DataType
    
    @abstractmethod
    def get_string(self) -> str:
        ...


@dataclass
class NFPythonToolProcessInput(NFProcessInput):

    def get_string(self) -> str:
        return f'path {self.name}'


@dataclass
class NFValProcessInput(NFProcessInput):

    def get_string(self) -> str:
        return f'val {self.name}'


@dataclass
class NFPathProcessInput(NFProcessInput):
    presents_as: Optional[str]=None

    @property
    def stage_as(self) -> str:
        # presents_as takes precent
        if self.presents_as:
            outstr = f", stageAs: '{self.presents_as}'"
        
        # optional files need to be staged in subdir to avoid name clashes
        elif self.dtype.optional:
            outstr = f", stageAs: '{self.name}/*'"
        
        # files should be staged according to original naming
        else:
            outstr = f", stageAs: '*'"
        return outstr
    
    def get_string(self) -> str:
        return f'path {self.name}{self.stage_as}'


@dataclass
class NFTupleProcessInput(NFProcessInput):
    qualifiers: list[str]
    subnames: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for qual, subname in zip(self.qualifiers, self.subnames):
            out += f'{qual}({subname}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out

    def get_string(self) -> str:
        return f'tuple {self.fields}'