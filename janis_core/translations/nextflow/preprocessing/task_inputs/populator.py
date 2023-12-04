

import os
from typing import Any
from abc import ABC, abstractmethod

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core import WorkflowBuilder, Tool, CodeTool, PythonTool, CommandToolBuilder, TInput
from janis_core.types import File
from janis_core import settings

from ....common import TaskInputCollector
from ... import params
from ... import naming
from ... import task_inputs 
from .common import get_true_workflow_inputs



### COMPOSITION HELPER METHODS ###

def populate_code_file(tool: PythonTool) -> None:
    if not isinstance(tool, PythonTool):
        return 
    
    path = f'{settings.translate.nextflow.TEMPLATES_OUTDIR}{os.sep}{tool.id()}.py'
    params.add(
        task_id=tool.id(),
        tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE,
        subtype='sub_tool',
        janis_dtype=File(),
        default=path
    )

    # update task inputs
    task_inputs.update(
        tool_id=tool.id(), 
        dstype_str='task_input', 
        tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE, 
        value=settings.translate.nextflow.PYTHON_CODE_FILE
    )

def populate_scripts(tool: CommandToolBuilder) -> None:
    # for CommandTools, need to have a task input for each script in files_to_create.
    # *unless translation is from Galaxy, where these will already be inputs. 
    if settings.ingest.SOURCE == 'galaxy' or not tool._files_to_create:
        return 
    
    for filename in tool._files_to_create.keys():  # type: ignore
        # get the file path to where the script will appear in the translation
        assert(isinstance(filename, str))
        path = os.path.join(settings.translate.nextflow.TEMPLATES_OUTDIR, filename)
        
        # generate a name for this input
        if len(tool._files_to_create) == 1:        # type: ignore
            name = 'script'
        else:
            name = naming.process.files_to_create_script(filename)

        # create param for nextflow.config & so we can get the param for process calls
        params.add(
            task_id=tool.id(),
            tinput_id=name,
            subtype='sub_tool',
            janis_dtype=File(),
            default=path
        )

        # update task inputs
        task_inputs.update(
            tool_id=tool.id(), 
            dstype_str='task_input', 
            tinput_id=name, 
            value=name
        )

def gen_task_input_value_process(tool: CommandToolBuilder | PythonTool, tinput_id: str) -> Any:
    tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
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

def gen_task_input_value_workflow(tool: WorkflowBuilder, tinput_id: str) -> Any:
    tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
    value = naming.process.generic(tinput)
    value = f'ch_{value}'
    return value


### TASK INPUT POPULATION CLASSES ###

class TaskInputsPopulator(ABC):
    
    @abstractmethod
    def populate(self) -> None:
        ...
    
    @abstractmethod
    def gen_task_input_value(self, tinput_id: str) -> Any:
        ...


class TaskInputsPopulatorToolIngest(TaskInputsPopulator):

    def __init__(self, tool: CommandToolBuilder | PythonTool) -> None:
        self.tool = tool
        if isinstance(tool, PythonTool):
            populate_code_file(tool)
        if isinstance(tool, CommandToolBuilder):
            populate_scripts(tool)

    def gen_task_input_value(self, tinput_id: str) -> Any:
        return gen_task_input_value_process(self.tool, tinput_id)
 
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
        ti_value = self.gen_task_input_value(tinput_id)
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    
    def update_as_static_input(self, tinput_id: str) -> None:
        tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
        ti_type = 'static' 
        ti_value = tinput.default
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    


class TaskInputsPopulatorCommandTool(TaskInputsPopulator):
    
    def __init__(self, tool: CommandToolBuilder, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
        self.tool = tool
        self.sources = sources
        self.main_wf = main_wf            
        
    def populate(self) -> None:
        populate_scripts(self.tool)

        if not task_inputs.exists(self.tool.id()):
            task_inputs.add_tool(self.tool.id())
        
        if settings.translate.MODE == 'extended':
            for tinput in self.tool.tool_inputs():
                self.update_as_task_input(tinput)

        elif settings.translate.MODE in ['skeleton', 'regular']:
            collector = TaskInputCollector(self.tool)
            collector.collect(self.main_wf)
            task_input_ids = [x for x in collector.histories.keys()]
            for tinput in self.tool.tool_inputs():
                if tinput.id() in task_input_ids:
                    self.update_as_task_input(tinput)
                elif tinput.default is not None:
                    self.update_as_static_input(tinput)
                else:
                    self.update_as_ignored_input(tinput)
        
    def update_as_task_input(self, tinput: TInput) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value(tinput.id())
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), ti_value)
    
    def update_as_static_input(self, tinput: TInput) -> None:
        ti_type = 'static'
        value = tinput.default
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), value)

    def update_as_ignored_input(self, tinput: TInput) -> None:
        ti_type = 'ignored'
        ti_value = None
        task_inputs.update(self.tool.id(), ti_type, tinput.id(), ti_value)
    
    def gen_task_input_value(self, tinput_id: str) -> Any:
        return gen_task_input_value_process(self.tool, tinput_id)



class TaskInputsPopulatorPythonTool(TaskInputsPopulator):
    
    def __init__(self, tool: PythonTool, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
        self.tool = tool
        self.sources = sources
        self.main_wf = main_wf

    def populate(self) -> None:
        populate_code_file(self.tool)

        if not task_inputs.exists(self.tool.id()):
            task_inputs.add_tool(self.tool.id())
        for tinput in self.tool.tool_inputs():
            self.update_as_task_input(tinput.id())

    def update_as_task_input(self, tinput_id: str) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value(tinput_id)
        task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    
    # def update_as_static_input(self, tinput_id: str) -> None:
    #     tinput = [x for x in self.tool.tool_inputs() if x.id() == tinput_id][0]
    #     ti_type = 'static' 
    #     ti_value = tinput.default
    #     task_inputs.update(self.tool.id(), ti_type, tinput_id, ti_value)
    
    def gen_task_input_value(self, tinput_id: str) -> Any:
        return gen_task_input_value_process(self.tool, tinput_id)



class TaskInputsPopulatorSubWorkflow(TaskInputsPopulator):
    
    def __init__(self, subwf: WorkflowBuilder, sources: dict[str, Any], main_wf: WorkflowBuilder) -> None:
        self.subwf = subwf
        self.sources = sources
        self.main_wf = main_wf

    def populate(self) -> None:
        if not task_inputs.exists(self.subwf.id()):
            task_inputs.add_tool(self.subwf.id())
        
        collector = TaskInputCollector(self.subwf)
        collector.collect(self.main_wf)
        task_input_ids = [x for x in collector.histories.keys()]

        for tinput in self.subwf.tool_inputs():
            if tinput.id() in task_input_ids:
                self.update_as_task_input(tinput.id())
            # elif tinput.default is not None:
            #     self.update_as_static_input(tinput)
            else:
                self.update_as_ignored_input(tinput.id())
            
    def update_as_task_input(self, tinput_id: str) -> None:
        ti_type = 'task_input'
        ti_value = self.gen_task_input_value(tinput_id)
        task_inputs.update(self.subwf.id(), ti_type, tinput_id, ti_value)

    def update_as_ignored_input(self, tinput_id: str) -> None:
        ti_type = 'ignored'
        ti_value = None
        task_inputs.update(self.subwf.id(), ti_type, tinput_id, ti_value)
    
    def gen_task_input_value(self, tinput_id: str) -> Any:
        return gen_task_input_value_workflow(self.subwf, tinput_id)



class TaskInputsPopulatorMainWorkflow(TaskInputsPopulator):
    
    def __init__(self, main_wf: WorkflowBuilder) -> None:
        self.main_wf = main_wf

    def populate(self) -> None:
        if settings.translate.MODE == 'extended':
            self.populate_task_inputs_extended()
        elif settings.translate.MODE in ['skeleton', 'regular']:
            self.populate_task_inputs_pruned()

    def populate_task_inputs_extended(self) -> None:
        all_tinput_ids = set([x.id() for x in self.main_wf.tool_inputs()])
        for tinput_id in all_tinput_ids:
            ti_type = 'param'
            tinput = [x for x in self.main_wf.tool_inputs() if x.id() == tinput_id][0]
            subtype = 'main_workflow'
            param = params.register(tinput, task_id=self.main_wf.id(), subtype=subtype)
            value = f'params.{param.name}' 
            task_inputs.update(self.main_wf.id(), ti_type, tinput_id, value)

    def populate_task_inputs_pruned(self) -> None:
        all_tinput_ids = set([x.id() for x in self.main_wf.tool_inputs()])
        param_tinput_ids = get_true_workflow_inputs(self.main_wf)
        ignored_tinput_ids = all_tinput_ids - param_tinput_ids

        # param inputs
        for tinput_id in param_tinput_ids:
            ti_type = 'param'
            tinput = [x for x in self.main_wf.tool_inputs() if x.id() == tinput_id][0]
            subtype = 'main_workflow'
            param = params.register(tinput, task_id=self.main_wf.id(), subtype=subtype)
            value = f'params.{param.name}'
            task_inputs.update(self.main_wf.id(), ti_type, tinput_id, value)
        
        # ignored inputs
        for tinput_id in ignored_tinput_ids:
            ti_type = 'ignored'
            value = None
            task_inputs.update(self.main_wf.id(), ti_type, tinput_id, value)
    
    def gen_task_input_value(self, tinput_id: str) -> Any:
        raise NotImplementedError



