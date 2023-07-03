

import os
from typing import Any
from abc import ABC, abstractmethod

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core import WorkflowBuilder, Tool, PythonTool, CommandTool, TInput
from janis_core.types import File
from janis_core import settings

from ... import params
from ... import naming
from ... import task_inputs 
from ....common import TaskInputCollector
from ....common import get_step_referenced_tinputs


class TaskInputsPopulator(ABC):
    
    def __init__(self, tool: Tool) -> None:
        self.tool = tool
        self.populate_code_file()
        self.populate_scripts()

    def populate_code_file(self) -> None:
        # pythontool gets extra code_file input before normal inputs
        if isinstance(self.tool, PythonTool):
            path = f'{settings.translate.nextflow.BASE_OUTDIR}{os.sep}{settings.translate.nextflow.TEMPLATES_OUTDIR}{os.sep}{self.tool.id()}.py'
            params.add(
                task_id=self.tool.id(),
                tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE,
                subtype='sub_tool',
                # name_override=self.tool.id(),
                janis_dtype=File(),
                default=path
            )

            # update task inputs
            task_inputs.update(
                tool_id=self.tool.id(), 
                dstype_str='task_input', 
                tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE, 
                value=settings.translate.nextflow.PYTHON_CODE_FILE
            )
    
    def populate_scripts(self) -> None:
        # for CommandTools, need to have a task input for each script in files_to_create.
        # *unless translation is from Galaxy, where these will already be inputs. 
        if isinstance(self.tool, CommandTool) and not settings.ingest.SOURCE == 'galaxy':
            if self.tool._files_to_create:                          # type: ignore
                for filename in self.tool._files_to_create.keys():  # type: ignore
                    
                    # get the file path to where the script will appear in the translation
                    assert(isinstance(filename, str))
                    path = os.path.join(settings.translate.nextflow.BASE_OUTDIR, settings.translate.nextflow.TEMPLATES_OUTDIR, filename)
                    
                    # generate a name for this input
                    if len(self.tool._files_to_create) == 1:        # type: ignore
                        name = 'script'
                    else:
                        name = naming.process.files_to_create_script(filename)

                    # create param for nextflow.config & so we can get the param for process calls
                    params.add(
                        task_id=self.tool.id(),
                        tinput_id=name,
                        subtype='sub_tool',
                        janis_dtype=File(),
                        default=path
                    )

                    # update task inputs
                    task_inputs.update(
                        tool_id=self.tool.id(), 
                        dstype_str='task_input', 
                        tinput_id=name, 
                        value=name
                    )
            
    @abstractmethod
    def populate(self) -> None:
        ...
        
    def gen_task_input_value_process(self, tinput_id: str) -> Any:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        # is_duplicate = self.duplicate_datatype_exists(tinput)
        dtt = utils.get_dtt(tinput.intype)

        if dtt == DTypeType.SECONDARY_ARRAY:
            value = naming.process.secondaries_array(tinput)
        elif dtt == DTypeType.SECONDARY:
            value = naming.process.generic(tinput)
        elif dtt == DTypeType.FILE_PAIR_ARRAY:
            value = naming.process.file_pair_array(tinput)
        elif dtt == DTypeType.FILE_PAIR:
            value = naming.process.file_pair(tinput)
        else:
            value = naming.process.generic(tinput)
        
        return value



class TaskInputsPopulatorToolMode(TaskInputsPopulator):

    def __init__(self, tool: Tool) -> None:
        super().__init__(tool)

    def populate(self) -> None:
        if not task_inputs.exists(self.tool.id()):
            task_inputs.add_tool(self.tool.id())
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

    def __init__(self, tool: Tool, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
        super().__init__(tool)
        self.sources = sources
        self.collector = TaskInputCollector(tool)
        self.collector.collect(main_wf)
    
    def populate(self) -> None:
        if not task_inputs.exists(self.tool.id()):
            task_inputs.add_tool(self.tool.id())
        if settings.translate.MODE in ['skeleton', 'regular']:
            task_input_ids = get_step_referenced_tinputs(self.collector)
        else:
            task_input_ids = [x.id() for x in self.tool.tool_inputs()]
        
        for tinput in self.tool.tool_inputs(): 
            if tinput.id() in task_input_ids:
                self.update_as_task_input(tinput)
            elif tinput.default is not None:
                self.update_as_static_input(tinput) 
            else:
                self.update_as_ignored_input(tinput) 

    ### helper methods
    def update_as_task_input(self, tinput: TInput) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value(tinput.id())
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), ti_value)
    
    def gen_task_input_value(self, tinput_id: str) -> None:
        if isinstance(self.tool, WorkflowBuilder):
            return self.gen_task_input_value_workflow(tinput_id)
        else:
            return self.gen_task_input_value_process(tinput_id)
    
    def gen_task_input_value_workflow(self, tinput_id: str) -> Any:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        value = naming.process.generic(tinput)
        value = f'ch_{value}'
        return value
    
    def update_as_static_input(self, tinput: TInput) -> None:
        ti_type = 'static'
        value = tinput.default
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), value)

    def update_as_ignored_input(self, tinput: TInput) -> None:
        ti_type = 'ignored'
        ti_value = None
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), ti_value)
