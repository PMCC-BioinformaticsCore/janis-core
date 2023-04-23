



import os
from typing import Any
from abc import ABC, abstractmethod

from janis_core import translation_utils as utils
from janis_core import Workflow, TInput, Tool, PythonTool
from janis_core.types import DataType, File
from janis_core import settings

from ... import params
from ... import naming
from ... import task_inputs

from .categories import TaskInputsCategoriser



class TaskInputsPopulator(ABC):
    
    def __init__(self, tool: Tool) -> None:
        self.tool = tool
        self.task_inputs: set[str] = set()
        self.param_inputs: set[str] = set()
        self.static_inputs: set[str] = set()
        self.ignored_inputs: set[str] = set()
        self.populate_code_file()

    def populate_code_file(self) -> None:
        # pythontool gets extra code_file input before normal inputs
        if isinstance(self.tool, PythonTool):
            path = f'{settings.translate.nextflow.BASE_OUTDIR}/{settings.translate.nextflow.TEMPLATES_OUTDIR}/{self.tool.id()}.py'
            # path = f'{os.getcwd()}/templates/{self.tool.id()}.py'
            # create param for nextflow.config & so we can get the param for process calls
            param = params.add(
                task_id=self.tool.id(),
                tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE_SYMBOL,
                subtype='sub_tool',
                # name_override=self.tool.id(),
                janis_dtype=File(),
                default=path
            )
            task_inputs.update(
                tool_id=self.tool.id(), 
                dstype_str='task_input', 
                tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE_SYMBOL, 
                value=settings.translate.nextflow.PYTHON_CODE_FILE_SYMBOL
            )

    @abstractmethod
    def populate(self) -> None:
        ...
        
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
            if tinput.id() in self.task_inputs:
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
               
    def gen_task_input_value_process(self, tinput_id: str) -> Any:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        dtype: DataType = tinput.intype  # type: ignore
        is_duplicate = self.duplicate_datatype_exists(tinput)
        
        if utils.is_array_secondary_type(dtype):
            value = naming.process.secondaries_array(tinput, duplicate_datatype_exists=is_duplicate)
        elif utils.is_secondary_type(dtype):
            value = naming.process.secondaries(tinput, duplicate_datatype_exists=is_duplicate)
        else:
            value = naming.process.generic(tinput)
        
        return value



class TaskInputsPopulatorToolMode(TaskInputsPopulator):

    def __init__(self, tool: Tool) -> None:
        super().__init__(tool)

    def populate(self) -> None:
        for tinput in self.tool.tool_inputs():
            if tinput.default is not None:
                self.update_as_static_input(tinput.id())
            else:
                self.update_as_task_input(tinput.id())

    def update_as_task_input(self, tinput_id: str) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value_process(tinput_id)
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    
    def update_as_static_input(self, tinput_id: str) -> None:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        ti_type = 'static'
        ti_value = tinput.default
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)



class TaskInputsPopulatorWorkflowMode(TaskInputsPopulator):

    def __init__(self, tool: Tool, sources: dict[str, Any], main_wf: Workflow) -> None:
        super().__init__(tool)
        self.sources = sources
        self.main_wf = main_wf
    
    def populate(self) -> None:
        self.update_categories()

        for tinput_id in self.task_inputs: 
            self.update_as_task_input(tinput_id)
        for tinput_id in self.param_inputs: 
            self.update_as_param_input(tinput_id)
        for tinput_id in self.static_inputs:
            try:
                self.update_as_static_input(tinput_id)
            except Exception:
                self.update_as_task_input(tinput_id)
        for tinput_id in self.ignored_inputs: 
            self.update_as_ignored_input(tinput_id)

    def update_categories(self) -> None:
        categoriser = TaskInputsCategoriser(self.tool, self.main_wf)
        categoriser.categorise()
        self.task_inputs = categoriser.task_inputs
        self.param_inputs = categoriser.param_inputs
        self.static_inputs = categoriser.static_inputs
        self.ignored_inputs = categoriser.ignored_inputs

    ### helper methods
    def update_as_task_input(self, tinput_id: str) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value(tinput_id)
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    
    def gen_task_input_value(self, tinput_id: str) -> None:
        if isinstance(self.tool, Workflow):
            return self.gen_task_input_value_workflow(tinput_id)
        else:
            return self.gen_task_input_value_process(tinput_id)
    
    def gen_task_input_value_workflow(self, tinput_id: str) -> Any:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        value = naming.process.generic(tinput)
        value = f'ch_{value}'
        return value

    def update_as_param_input(self, tinput_id: str) -> None:
        ti_type = 'param'
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        
        # task subtype
        if self.tool.id() == self.main_wf.id():
            subtype = 'main_workflow'
        elif isinstance(self.tool, Workflow):
            subtype = 'sub_workflow'
        else:
            subtype = 'sub_tool'
        
        param = params.register(tinput, self.tool.id(), subtype)
        value = f'params.{param.name}'
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
