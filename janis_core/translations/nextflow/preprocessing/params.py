

import os
from copy import deepcopy

from janis_core.workflow.workflow import Workflow, InputNode, StepNode
from janis_core.types import File
from janis_core import PythonTool
from janis_core import translation_utils as utils
from janis_core import settings

from .channels import get_channel_inputs_to_register
from .. import params
from ..scope import Scope
    

def register_params(wf: Workflow) -> None:
    """
    register param(s) for each workflow input. 
    channel(s) may also be registered if necessary.
    """
    scope = Scope()
    param_inputs = get_param_inputs_to_register(wf, scope)
    channel_inputs = get_channel_inputs_to_register(wf, scope)

    for inp in wf.input_nodes.values():
        is_channel_input = True if inp.id() in channel_inputs else False
        if inp.id() in param_inputs:
            ParamRegistrationHelper(inp, scope, is_channel_input).register()
    
    # repeat for nested workflows (subworkflows)
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        if isinstance(step.tool, PythonTool):
            register_params_python_tool(step, current_scope)


def register_params_python_tool(step: StepNode, current_scope: Scope) -> None:
    """A param will be registered for the code_file of each PythonTool."""
    default = _get_code_file_path(step.tool)
    params.add(
        janis_tag='code_file',
        scope=current_scope,
        default=default,
        is_channel_input=False,
        janis_dtype=File(),
        janis_uuid=None,
    )


# helper classes 

class ParamRegistrationHelper:
    def __init__(self, inp: InputNode, scope: Scope, is_channel_input: bool) -> None:
        self.inp = inp
        self.scope = scope
        self.is_channel_input = is_channel_input

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
        params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            is_channel_input=True,
            janis_dtype=self.inp.datatype,
            janis_uuid=self.inp.uuid
        )
    
    def register_param_secondaries(self) -> None:
        params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            is_channel_input=True,
            janis_dtype=self.inp.datatype,
            janis_uuid=self.inp.uuid
        )
    
    def register_param(self) -> None:
        params.add(
            janis_tag=self.inp.id(),
            scope=self.scope,
            default=self.inp.default if self.inp.default is not None else None,
            is_channel_input=self.is_channel_input,
            janis_dtype=self.inp.datatype,
            janis_uuid=self.inp.uuid,
        )





# helper functions 
# these identify which workflow inputs should have a param / channel created. 

def get_param_inputs_to_register(wf: Workflow, scope: Scope) -> set[str]:
    if scope.labels == [settings.translate.nextflow.NF_MAIN_NAME]:
        items: set[str] = {x.id() for x in wf.input_nodes.values()}
    else:
        items: set[str] = {x.id() for x in wf.input_nodes.values()} - get_channel_inputs_to_register(wf, scope)
    return items

def _get_code_file_path(tool: PythonTool) -> str:
    basedir = settings.translate.nextflow.BASE_OUTDIR
    subfolder = settings.translate.nextflow.CODE_FILES_OUTDIR
    filename = tool.id()
    filepath = os.path.join(basedir, subfolder, filename)
    filepath += '.py'
    return filepath


