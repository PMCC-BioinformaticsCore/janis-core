

from typing import Any

from janis_core import TInput
from janis_core.types import DataType
from janis_core import translation_utils as utils

from .main import add
from .main import Param


def register(tinput: TInput, task_id: str, subtype: str) -> Param:
    manager = ParamRegistrationManager(tinput, task_id, subtype)
    return manager.register()


class ParamRegistrationManager:
    
    def __init__(self, tinput: TInput, task_id: str, subtype: str) -> None:
        self.tinput = tinput
        self.task_id = task_id
        self.subtype = subtype
    
    def register(self) -> Param:
        """registers param for each wf input which requires a param."""
        new_param = add(
            task_id=self.task_id,
            tinput_id=self.tinput.id(),
            subtype=self.subtype,
            janis_dtype=self.tinput.intype,
            default=self.default,
        )
        return new_param

    @property
    def tformat(self) -> str:
        # secondaries array
        if utils.is_array_secondary_type(self.tinput.intype):
            return 'secondary_array'
        # secondaries
        elif utils.is_secondary_type(self.tinput.intype):
            return 'secondary'
        # file_pair array
        elif utils.is_array_file_pair_type(self.tinput.intype):
            return 'file_pair_array'
        # file_pair
        elif utils.is_file_pair_type(self.tinput.intype):
            return 'file_pair'
        # TODO array file pairs?
        # file array
        elif self.tinput.intype.is_array():
            return 'generic_array'
        # anything else
        else:
            return 'generic'
        
    @property 
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.tinput.intype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property 
    def default(self) -> Any:

        # default given
        if self.tinput.default is not None:
            return self.tinput.default
        
        # no default, but optional file type
        if utils.is_file_type(self.tinput.intype) and self.tinput.intype.optional:
            if self.tformat == 'secondary_array':
                return self.get_file_default_secondary_array()
            
            elif self.tformat == 'secondary':
                return self.get_file_default_secondary()
            
            elif self.tformat == 'file_pair_array':
                return self.get_file_default_file_pair_array()
            
            elif self.tformat == 'file_pair':
                return self.get_file_default_file_pair()
            
            elif self.tformat == 'generic_array':
                return self.get_file_default_generic_array()
            
            elif self.tformat == 'generic':
                return self.get_file_default_generic()
            
            else:
                raise RuntimeError
        
        return None
    
    def get_file_default_secondary_array(self) -> list[list[str]]:
        basetype = utils.get_base_type(self.tinput.intype)
        num_files = len(utils.get_extensions(basetype))
        internal_structure = ['NO_FILE'] * num_files
        default = [internal_structure]
        return default
    
    def get_file_default_secondary(self) -> list[str]:
        basetype = utils.get_base_type(self.tinput.intype)
        num_files = len(utils.get_extensions(basetype))
        default = ['NO_FILE'] * num_files
        return default
    
    def get_file_default_file_pair_array(self) -> list[list[str]]:
        default = [['NO_FILE', 'NO_FILE']]
        return default
    
    def get_file_default_file_pair(self) -> list[str]:
        default = ['NO_FILE', 'NO_FILE']
        return default
    
    def get_file_default_generic_array(self) -> list[str]:
        default = ['NO_FILE']
        return default
    
    def get_file_default_generic(self) -> str:
        default = 'NO_FILE'
        return default
            
