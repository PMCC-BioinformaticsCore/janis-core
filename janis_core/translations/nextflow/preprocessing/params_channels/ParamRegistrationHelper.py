

from janis_core.workflow.workflow import InputNode
from janis_core import translation_utils as utils

from ... import params
from ...scope import Scope
    


class ParamRegistrationHelper:
    def __init__(self, inp: InputNode, scope: Scope) -> None:
        self.inp = inp
        self.scope = scope

    def register(self) -> None:
        """registers param for each wf input which requires a param."""
        
        # secondaries array
        if utils.is_array_secondary_type(self.inp.datatype):
            self.register_param_secondaries_array()
        
        # secondaries
        elif utils.is_secondary_type(self.inp.datatype):
            self.register_param_secondaries()

        # anything else
        else:
            self.register_param()
    
    def register_param_secondaries_array(self) -> None:  
        # @secondariesarray
        new_param = params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            janis_dtype=self.inp.datatype,
        )
        params.create_link(janis_uuid=self.inp.uuid, param=new_param)
    
    def register_param_secondaries(self) -> None:
        new_param = params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            janis_dtype=self.inp.datatype,
        ) 
        params.create_link(janis_uuid=self.inp.uuid, param=new_param)
    
    def register_param(self) -> None:
        new_param = params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            default=self.inp.default if self.inp.default is not None else None,
            janis_dtype=self.inp.datatype,
        )
        params.create_link(janis_uuid=self.inp.uuid, param=new_param)


