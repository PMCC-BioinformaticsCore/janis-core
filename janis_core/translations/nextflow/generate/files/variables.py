
from typing import Optional
from enum import Enum, auto

from janis_core import settings
from janis_core import Workflow, TInput, DataType
from janis_core import translation_utils as utils

from ... import naming 
from ... import task_inputs 
from ...task_inputs import TaskInputType
from ...task_inputs import TaskInput
from ...model.files import NFVariableDefinition
from ...model.files import NFVariableDefinitionBlock
from ...model.workflow import NFWorkflow
from ...model.workflow import NFMainWorkflow


INDENT = settings.translate.nextflow.NF_INDENT
CROSS_CHANNEL_NAME = 'ch_cartesian_cross'






def gen_variables_block(nf_workflow: NFWorkflow, wf: Workflow) -> Optional[NFVariableDefinitionBlock]:
    # for each param
    # if it qualifies to create a file, create file variable declaration
    # return block if 1+ variable declarations 
    var_definitions: list[NFVariableDefinition] = []
    var_block: Optional[NFVariableDefinitionBlock] = None
    
    if isinstance(nf_workflow, NFMainWorkflow):
        for tinput in wf.tool_inputs():
            if task_inputs.exists(wf.id(), tinput):
                task_input = task_inputs.get(wf.id(), tinput)
                if task_input.ti_type in (TaskInputType.PARAM, TaskInputType.TASK_INPUT):
                    # only optional file types get declared as file() variables
                    if tinput.intype.optional and utils.is_file_type(tinput.intype):
                        generator = VariableDefinitionGenerator(tinput, task_input)
                        var_def = generator.generate()
                        var_definitions.append(var_def)
                        
    if var_definitions:
        var_block = NFVariableDefinitionBlock(var_definitions)
    
    return var_block



class VariableValueType(Enum):
    FILE_OPT            = auto()
    FILEPAIR_OPT        = auto()
    SECONDARY_OPT       = auto()
    FILE_ARR_OPT        = auto()
    FILEPAIR_ARR_OPT    = auto()
    SECONDARY_ARR_OPT   = auto()


class VariableDefinitionGenerator:
    def __init__(self, tinput: TInput, task_input: TaskInput) -> None:
        self.tinput = tinput
        self.task_input = task_input

    @property
    def dtype(self) -> DataType:
        return self.tinput.intype  # type: ignore
    
    @property
    def varname(self) -> str:
        return naming.constructs.gen_varname_file(self.tinput.id())
    
    @property
    def paramname(self) -> str:
        assert(isinstance(self.task_input.value, str))
        return self.task_input.value
    
    @property
    def vartype(self) -> VariableValueType:
        if utils.is_array_secondary_type(self.dtype):
            return VariableValueType.SECONDARY_ARR_OPT
        elif utils.is_secondary_type(self.dtype):
            return VariableValueType.SECONDARY_OPT
        elif utils.is_array_file_pair_type(self.dtype):
            return VariableValueType.FILEPAIR_ARR_OPT
        elif utils.is_file_pair_type(self.dtype):
            return VariableValueType.FILEPAIR_OPT
        elif self.dtype.is_array():
            return VariableValueType.FILE_ARR_OPT
        else:
            return VariableValueType.FILE_OPT
    
    @property
    def varvalue(self) -> str:
        if self.vartype == VariableValueType.FILE_OPT:
            value =  f'file( {self.paramname} )'
        
        elif self.vartype == VariableValueType.FILEPAIR_OPT:
            value = f'{self.paramname}.each {{ value -> file(value) }}'

        elif self.vartype == VariableValueType.SECONDARY_OPT:
            value = f'{self.paramname}.each {{ value -> file(value) }}'

        elif self.vartype == VariableValueType.FILE_ARR_OPT:
            value = f'{self.paramname}.each {{ value -> file(value) }}'

        elif self.vartype == VariableValueType.FILEPAIR_ARR_OPT:
            value = f'{self.paramname}.each {{ outer -> outer.each {{ inner -> file(inner) }} }}'

        elif self.vartype == VariableValueType.SECONDARY_ARR_OPT:
            value = f'{self.paramname}.each {{ outer -> outer.each {{ inner -> file(inner) }} }}'

        return value

    def generate(self) -> NFVariableDefinition:
        return NFVariableDefinition(
            name=self.varname,
            value=self.varvalue, 
            dtype=self.dtype
        )
        
