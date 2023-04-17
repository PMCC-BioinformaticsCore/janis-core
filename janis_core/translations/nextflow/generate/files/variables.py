
from typing import Optional


from janis_core import settings
from janis_core import Workflow, TInput
from janis_core import translation_utils as utils

from ... import naming 
from ... import task_inputs 
from ...task_inputs import TaskInputType
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
        for input_node in wf.tool_inputs():
            if task_inputs.exists(wf.id(), input_node):
                task_input = task_inputs.get(wf.id(), input_node)
                if task_input.ti_type in (TaskInputType.PARAM, TaskInputType.TASK_INPUT):
                    if input_node.intype.optional:
                        f_name = naming.constructs.gen_varname_file(input_node.id())
                        f_source = _get_source(input_node, wf)
                        f_def = NFVariableDefinition(f_name, f_source, input_node.intype)
                        var_definitions.append(f_def)

    if var_definitions:
        var_block = NFVariableDefinitionBlock(var_definitions)
    
    return var_block


def _get_source(tinput: TInput, wf: Workflow) -> str:
    task_input = task_inputs.get(wf.id(), tinput)
    param_name = task_input.value
    # @secondaryarrays
    if utils.is_array_secondary_type(tinput.intype):
        src = f'{param_name}.flatten()'
    else:
        src = f'{param_name}'
    return src
