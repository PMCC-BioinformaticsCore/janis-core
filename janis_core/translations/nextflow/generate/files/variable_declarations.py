
from typing import Optional

from janis_core import settings
from janis_core import Workflow, TInput, DataType
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ... import naming 
from ... import task_inputs 
from ...task_inputs import TaskInputType
from ...task_inputs import TaskInput
from ...model.files import NFVariableDefinition
from ...model.files import NFVariableDefinitionBlock
from ...model.workflow import NFWorkflow
from ...model.workflow import NFMainWorkflow
from ...nulls import add_file_cast

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
            if should_create_variable_definition(tinput, wf):
                task_input = task_inputs.get(wf.id(), tinput)
                generator = VariableDefinitionGenerator(tinput, task_input)
                var_def = generator.generate()
                var_definitions.append(var_def)
                        
    if var_definitions:
        var_block = NFVariableDefinitionBlock(var_definitions)
    
    return var_block


def should_create_variable_definition(tinput: TInput, wf: Workflow) -> bool:
    if not task_inputs.exists(wf.id(), tinput):
        return False
    
    task_input = task_inputs.get(wf.id(), tinput)
    if task_input.ti_type in (TaskInputType.STATIC, TaskInputType.IGNORED, TaskInputType.LOCAL):
        return False
    
    dtt = utils.get_dtt(tinput.intype)
    if dtt not in [
        DTypeType.SECONDARY_ARRAY,
        DTypeType.SECONDARY,
        DTypeType.FILE_PAIR_ARRAY,
        DTypeType.FILE_PAIR,
        DTypeType.FILE_ARRAY,
        DTypeType.FILE,
    ]:
        return False

    if not tinput.intype.optional:
        return False
    
    return True


class VariableDefinitionGenerator:
    def __init__(self, tinput: TInput, task_input: TaskInput) -> None:
        self.tinput = tinput
        self.task_input = task_input

    @property
    def dtype(self) -> DataType:
        return self.tinput.intype  # type: ignore
    
    @property
    def var_name(self) -> str:
        return naming.constructs.gen_varname_file(self.tinput.id())
    
    @property
    def param_name(self) -> str:
        assert(isinstance(self.task_input.value, str))
        return self.task_input.value
    
    @property
    def value(self) -> str:
        dtt = utils.get_dtt(self.dtype)
        if dtt not in [
            DTypeType.SECONDARY_ARRAY,
            DTypeType.SECONDARY,
            DTypeType.FILE_PAIR_ARRAY,
            DTypeType.FILE_PAIR,
            DTypeType.FILE_ARRAY,
            DTypeType.FILE,
        ]:
            raise RuntimeError(f'Unsupported file type: {self.dtype}')

        value = self.param_name
        value = add_file_cast(self.dtype, value)
        return value

    def generate(self) -> NFVariableDefinition:
        return NFVariableDefinition(
            name=self.var_name,
            value=self.value, 
            dtype=self.dtype
        )
        
