


from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import DataType
# from janis_core import translation_utils as utils


@dataclass
class NFProcessInput(ABC):
    name: str
    tinput_id: str
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
        if self.presents_as:
            expr = f", stageAs: '{self.presents_as}'"

        # # optional secondary arrays
        # elif utils.is_secondary_array_type(self.dtype) and self.dtype.optional:
        #     expr = f", stageAs: '{self.name}??/*'"

        # # optional file pair arrays
        # elif utils.is_file_pair_array_type(self.dtype) and self.dtype.optional:
        #     expr = f", stageAs: '{self.name}??/*'"

        # # optional file arrays
        # elif utils.is_file_array_type(self.dtype) and self.dtype.optional:
        #     # expr = f", stageAs: '{self.name}??/*'"
        #     expr = ''
        
        # # optional files
        # elif utils.is_file_type(self.dtype) and self.dtype.optional:
        #     expr = f", stageAs: '{self.name}/*'"
        
        else:
            expr = ''
            # raise RuntimeError('CHECK ME: NFPathProcessInput')

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