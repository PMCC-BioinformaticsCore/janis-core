
from typing import Tuple

from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowInput
from ..TextRender import TextRender
#from ..formatting import format_typestr

def note_snippet(commenter: str) -> str:
    return f"""{commenter} WORKFLOW INPUTS
{commenter} This file contains workflow inputs which need to be provided by the user.
{commenter} The names of workflow inputs (below) will appear in workflow.py where they are used. 

{commenter} NULL VALUES
{commenter} If a value should be left null, ensure the tool input is optional. 
{commenter} The engine will throw an error otherwise.  
"""

# def note_snippet(commenter: str) -> str:
#     return f"""{commenter} This file contains workflow inputs which need to be provided by the user.
# {commenter} Organised as follows: 
# {commenter}     1. input data for the workflow
# {commenter}     2. runtime values for each step

# {commenter} Null values must be replaced by the user to run the workflow. 
# {commenter} VIGNETTES
# """

class InputsText(TextRender):
    def __init__(self, entity: Workflow, file_format: str='yaml'):
        super().__init__()
        self.entity = entity
        self.file_format = file_format

    @property
    def imports(self) -> list[Tuple[str, str]]:
        raise NotImplementedError()

    def render(self) -> str:
        if self.file_format == 'yaml':
            return to_yaml(self.entity.inputs)
        else:
            raise RuntimeError('wrong format')



YAML_NONE = 'null'
YAML_COMMENTER = '#'

def to_yaml(inputs: list[WorkflowInput]) -> str:
    out_str: str = ''
    out_str += '\n'
    out_str += f'{note_snippet(commenter=YAML_COMMENTER)}'
    out_str += '\n'

    for winp in inputs:
        if winp.value is None:
            winp.value = YAML_NONE
        #out_str += f'{winp.tag}: {winp.value}  {YAML_COMMENTER} [{format_typestr(winp)}]\n'
        out_str += f'{winp.tag}: {winp.value}\n'

    return out_str


# future formats here
