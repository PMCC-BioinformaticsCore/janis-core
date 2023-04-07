



from typing import Any

from janis_core import translation_utils as utils
from janis_core import Workflow, TInput, Tool
from janis_core.types import DataType

from ... import naming
from ... import task_inputs

from .categories import TaskInputsCategoriser



class TaskInputsPopulator:
    def __init__(self, tool: Tool, sources: dict[str, Any], main_wf: Workflow) -> None:
        self.tool = tool
        self.sources = sources
        self.main_wf = main_wf
        self.categoriser = TaskInputsCategoriser(self.tool, self.main_wf)
        self.categoriser.categorise()
        print()
    
    def populate(self) -> None:
        for tinput_id in self.categoriser.task_inputs: 
            self.update_as_task_input(tinput_id)
        for tinput_id in self.categoriser.param_inputs: 
            self.update_as_param_input(tinput_id)
        for tinput_id in self.categoriser.static_inputs: 
            self.update_as_static_input(tinput_id)
        for tinput_id in self.categoriser.ignored_inputs: 
            self.update_as_ignored_input(tinput_id)

    ### helper methods
    def update_as_task_input(self, tinput_id: str) -> None:
        ti_type = 'task_input'
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        dtype: DataType = tinput.intype  # type: ignore
        is_duplicate = self.duplicate_datatype_exists(tinput)
        
        if utils.is_array_secondary_type(dtype):
            value = naming.process.secondaries_array(tinput, duplicate_datatype_exists=is_duplicate)
        elif utils.is_secondary_type(dtype):
            value = naming.process.secondaries(tinput, duplicate_datatype_exists=is_duplicate)
        else:
            value = naming.process.generic(tinput)  
        task_inputs.update(self.tool.id(), ti_type, tinput_id, value)
    
    def update_as_param_input(self, tinput_id: str) -> None:
        ti_type = 'param'
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        pname = naming.constructs.gen_varname_param(
            tinput_id=tinput_id, 
            task_id=self.tool.id(), 
            name_override=None,
            dtype=tinput.intype
        )
        value = f'params.{pname}'
        task_inputs.update(self.tool.id(), ti_type, tinput_id, value)
    
    def update_as_static_input(self, tinput_id: str) -> None:
        ti_type = 'static'
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        src = self.sources[tinput.id()]
        node = utils.resolve_node(src)
        value = node.default
        task_inputs.update(self.tool.id(), ti_type, tinput_id, value)
    
    def update_as_ignored_input(self, tinput_id: str) -> None:
        ti_type = 'ignored'
        value = None
        task_inputs.update(self.tool.id(), ti_type, tinput_id, value)

    def duplicate_datatype_exists(self, inp: TInput) -> bool:
        """
        check if another TInput has the same dtype as this TInput.
        only checking for secondary and array secondary datatype duplicates.
        """
        rtype = inp.intype  # type: ignore
        rbasetype = utils.get_base_type(rtype)  # type: ignore
        
        for tinput in self.tool.tool_inputs():
            # dont check tinput against itself
            if tinput.id() == inp.id():
                continue
            
            # only check against other task inputs
            if tinput.id() in self.categoriser.task_inputs:
                # check if types match
                qtype = tinput.intype  # type: ignore
                qbasetype = utils.get_base_type(qtype)  # type: ignore

                if utils.is_array_secondary_type(rtype) and utils.is_array_secondary_type(qtype):  # type: ignore
                    if type(rbasetype) == type(qbasetype):
                        return True
                
                elif utils.is_secondary_type(rtype) and utils.is_secondary_type(qtype):  # type: ignore
                    if type(rtype) == type(qtype):  # type: ignore
                        return True
        return False
               
