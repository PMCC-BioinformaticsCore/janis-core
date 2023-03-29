

from typing import Any

from janis_core.workflow.workflow import InputNode
from janis_core.types import DataType, File, Filename, Directory
from janis_core import translation_utils as utils

from ... import params
from ...scope import Scope
    

no_file_count: int = 0


class ParamRegistrationHelper:
    def __init__(self, inp: InputNode, scope: Scope) -> None:
        self.inp = inp
        self.scope = scope
    
    def register(self) -> None:
        """registers param for each wf input which requires a param."""
        new_param = params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            default=self.default,
            janis_dtype=self.inp.datatype,
        )
        params.create_link(janis_uuid=self.inp.uuid, param=new_param)

    @property
    def tformat(self) -> str:
        # secondaries array
        if utils.is_array_secondary_type(self.inp.datatype):
            return 'secondary_array'
        # secondaries
        elif utils.is_secondary_type(self.inp.datatype):
            return 'secondary'
        # file_pair
        elif utils.is_file_pair_type(self.inp.datatype):
            return 'file_pair'
        # TODO array file pairs
        # file array
        elif self.inp.datatype.is_array():
            return 'generic_array'
        # anything else
        else:
            return 'generic'

    @property 
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.inp.datatype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property 
    def default(self) -> Any:

        # default given
        if self.inp.default is not None:
            return self.inp.default
        
        # no default, but optional type
        if isinstance(self.basetype, (File, Filename, Directory)) and self.inp.datatype.optional:
            if self.tformat == 'secondary_array':
                return self.get_default_secondary_array()
            
            elif self.tformat == 'secondary':
                return self.get_default_secondary()
            
            elif self.tformat == 'file_pair':
                return self.get_default_file_pair()
            
            elif self.tformat == 'generic_array':
                return self.get_default_generic_array()
            
            elif self.tformat == 'generic':
                return self.get_default_generic()
            
            else:
                raise RuntimeError
        
        return None
    
    def get_default_secondary_array(self) -> list[list[str]]:
        global no_file_count
        out: list[str] = []
        
        basetype = utils.get_base_type(self.inp.datatype)
        num_files = len(utils.get_extensions(basetype))
        for i in range(num_files):
            no_file_count += 1
            filename = f'NO_FILE{no_file_count}'
            out.append(filename)

        return [out]
    
    def get_default_secondary(self) -> list[str]:
        global no_file_count
        out: list[str] = []
        
        basetype = utils.get_base_type(self.inp.datatype)
        num_files = len(utils.get_extensions(basetype))
        for i in range(num_files):
            no_file_count += 1
            filename = f'NO_FILE{no_file_count}'
            out.append(filename)

        return out
    
    def get_default_file_pair(self) -> list[str]:
        global no_file_count
        out: list[str] = []
        for i in range(2):
            no_file_count += 1
            filename = f'NO_FILE{no_file_count}'
            out.append(filename)
        return out
    
    def get_default_generic_array(self) -> list[str]:
        global no_file_count
        no_file_count += 1
        out = [f'NO_FILE{no_file_count}']
        return out
    
    def get_default_generic(self) -> str:
        global no_file_count
        no_file_count += 1
        return f'NO_FILE{no_file_count}'
            

    
    # def register_param_secondaries_array(self) -> None:
    #     # @secondariesarray
    #     new_param = params.add(
    #         janis_tag=self.inp.id(),
    #         scope=self.scope,
    #         default=self.default,
    #         janis_dtype=self.inp.datatype,
    #     )
    #     params.create_link(janis_uuid=self.inp.uuid, param=new_param)
    
    # def register_param_secondaries(self) -> None:
    #     new_param = params.add(
    #         janis_tag=self.inp.id(),
    #         scope=self.scope,
    #         default=self.default,
    #         janis_dtype=self.inp.datatype,
    #     ) 
    #     params.create_link(janis_uuid=self.inp.uuid, param=new_param)
    
    # def register_param(self) -> None:
    #     new_param = params.add(
    #         janis_tag=self.inp.id(),
    #         scope=self.scope,
    #         default=self.default,
    #         janis_dtype=self.inp.datatype,
    #     )
    #     params.create_link(janis_uuid=self.inp.uuid, param=new_param)


