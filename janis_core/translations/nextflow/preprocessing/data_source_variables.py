


from typing import Optional, Any

from janis_core import (
    DataType,
    TInput,
)

from janis_core import translation_utils as utils
from janis_core import settings

from .. import params
from ..scope import Scope
from ..casefmt import to_case


from copy import deepcopy
from typing import Any

import janis_core.translation_utils as utils
from janis_core.workflow.workflow import Workflow
from janis_core import CommandTool, PythonTool, TInput
from janis_core import settings

from ..scope import Scope
from .. import naming
from .. import params
from ..process import data_sources



def register_ds_variables(entity: Workflow | CommandTool | PythonTool) -> None:
    """
    MAIN ENTRY POINT for this preprocessing task.

    for each CommandTool / PythonTool, decides the variable_name (a process_input, a param_input, or None ) 
    which feeds data for the tinput.
    the variable_name is how the particular tinput_id will be referenced inside its process. 
    """
    
    scope = Scope()

    if isinstance(entity, Workflow):
        register_ds_variables_workflow(scope, entity)
    elif isinstance(entity, (CommandTool | PythonTool)):  # type: ignore
        sources: dict[str, Any] = {}
        register_ds_variables_tool(scope, entity, sources)
    else:
        raise RuntimeError

def register_ds_variables_workflow(scope: Scope, wf: Workflow) -> None:
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        
        if isinstance(step.tool, CommandTool) or isinstance(step.tool, PythonTool):
            register_ds_variables_tool(current_scope, step.tool, step.sources)
        elif isinstance(step.tool, Workflow):
            register_ds_variables_workflow(current_scope, step.tool)
        else:
            raise NotImplementedError

def register_ds_variables_tool(scope: Scope, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> None:
    generator = ProcessDSVariableGenerator(scope, tool, sources)
    generator.generate()
    for tinput_id, variable_name in generator.names_map.items():
        data_sources.update_variables(scope, tinput_id, variable_name)


class ProcessDSVariableGenerator:
    """
    for each CommandTool / PythonTool, determines the internal variable_name for each tinput.id()
    the variable_name is how the particular tinput_id will be referenced inside the process. 
    """
    def __init__(self, scope: Scope, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> None:
        self.scope = scope
        self.tool = tool
        self.sources = sources

        # want to calculate this
        self.names_map: dict[str, Optional[str | list[str]]] = {}

    # helper properties
    @property
    def tinputs(self) -> list[TInput]:
        if isinstance(self.tool, CommandTool):
            return self.tool.tool_inputs()
        elif isinstance(self.tool, PythonTool):  # type: ignore
            return self.tool.inputs()
        else:
            raise NotImplementedError
    
    # helper methods
    def duplicate_datatype_exists(self, inp: TInput) -> bool:
        """
        check if another TInput has the same dtype as this TInput.
        only checking for secondary and array secondary datatype duplicates.
        """
        rtype = inp.intype  # type: ignore
        rbasetype = utils.get_base_type(rtype)  # type: ignore
        
        for tinput in self.tinputs:
            # dont check tinput against itself
            if tinput.id() == inp.id():
                continue
            
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
        
    def generate(self) -> None:
        for tinput in self.tinputs:
            name = self.gen_varname(tinput)
            self.names_map[tinput.id()] = name

    # this should be kept in ProcessVariableManager or something
    def gen_varname(self, inp: TInput) -> Optional[str | list[str]]:
        if inp.id() in data_sources.process_inputs(self.scope):
            return self.gen_varname_process_input(inp)
        elif inp.id() in data_sources.param_inputs(self.scope):
            return self.gen_varname_param_input(inp)
        elif inp.id() in data_sources.internal_inputs(self.scope):
            return None
        else:
            raise RuntimeError
        
    def gen_varname_process_input(self, inp: TInput) -> str | list[str]:
        # data fed via process input
        dtype: DataType = inp.intype  # type: ignore
        is_duplicate = self.duplicate_datatype_exists(inp)
        
        if utils.is_array_secondary_type(dtype):
            name = naming.process.secondaries_array(inp, duplicate_datatype_exists=is_duplicate)
        elif utils.is_secondary_type(dtype):
            name = naming.process.secondaries(inp, duplicate_datatype_exists=is_duplicate)
        else:
            name = naming.process.generic(inp)  
        return name
        
    def gen_varname_param_input(self, inp: TInput) -> str: 
        # data fed via global param
        src = self.sources[inp.id()]
        sel = src.source_map[0].source
        param = params.get(sel.input_node.uuid)
        pname = to_case(param.name, case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'params.{pname}'


