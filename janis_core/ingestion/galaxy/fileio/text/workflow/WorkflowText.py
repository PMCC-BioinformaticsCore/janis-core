



from __future__ import annotations
from typing import TYPE_CHECKING, Tuple

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.model.workflow import Workflow

from datetime import datetime
from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT

from .. import formatting
from .. import ordering

from ..TextRender import TextRender
from .WorkflowInputText import WorkflowInputText
from .WorkflowOutputText import WorkflowOutputText
from .StepText import StepText


### HELPER METHODS ### 

def note_snippet(workflow: Workflow) -> str:
    return f"""# NOTE
# This is an automated translation of the '{workflow.metadata.name}' version '{workflow.metadata.version}' workflow. 
# Translation was performed by the gxtool2janis program (in-development)
"""

def inputs_yaml_snippet() -> str:
    return f"""
\"\"\"
NOTE
Please provide values for these workflow inputs 
in inputs.yaml
\"\"\"
"""

def syspath_snippet() -> str: # ???
    return """import sys
sys.path.append('/home/grace/work/pp/gxtool2janis')
"""

def metadata_snippet(workflow: Workflow) -> str:
    contributors = ['gxtool2janis']
    return f"""metadata = WorkflowMetadata(
    short_documentation="{workflow.metadata.annotation}",
    contributors={contributors},
    keywords={workflow.metadata.tags},
    dateCreated="{datetime.today().strftime(JANIS_DATE_FMT)}",
    dateUpdated="{datetime.today().strftime(JANIS_DATE_FMT)}",
    version={workflow.metadata.version}
)
"""

def builder_snippet(workflow: Workflow) -> str:
    out_str: str = ''
    out_str += 'w = WorkflowBuilder(\n'
    out_str += f'\t"{workflow.tag}",\n'
    out_str += f'\tversion="{workflow.metadata.version}",\n'
    out_str += f'\tdoc="{workflow.metadata.annotation}",\n'
    out_str += f'\tmetadata=metadata\n'
    out_str += ')\n'
    return out_str

def translate_snippet() -> str:
    return f"""if __name__ == "__main__":
    w.translate("cwl", **args)
"""

core_imports = [
    ("janis_core", "WorkflowMetadata"),
    ("janis_core", "WorkflowBuilder"),
]

def label(text: str) -> str: 
    title = f'# {text}'
    border = f'# {"-" * (len(title) - 2)}'
    return f"""\
{border}
{title}
{border}"""

### MAIN CLASS ### 

class WorkflowText(TextRender):
    def __init__(
        self, 
        entity: Workflow, 
        render_note: bool=True,
        render_syspath: bool=True,
        render_translate: bool=True
    ):
        super().__init__()
        self.entity = entity
        self.render_note = render_note
        self.render_syspath = render_syspath
        self.render_translate = render_translate

    @property
    def imports(self) -> list[Tuple[str, str]]:
        imports: list[Tuple[str, str]] = []
        imports += core_imports
        # imports for each workflow input (eg datatypes)
        for winp in self.entity.inputs:
            imports += WorkflowInputText(winp).imports
        # imports for each translated tool
        for step in self.entity.steps:
            tool_id = step.metadata.wrapper.tool_id
            relative_path = f'tools.{tool_id}'
            tool_tag = step.tool.tag
            imports.append((relative_path, tool_tag))
        # imports used in steps
        for step in self.entity.steps:
            imports += StepText(step, self.entity).imports
        # imports for each output (mostly datatypes)
        for wout in self.entity.outputs:
            imports += WorkflowOutputText(wout).imports
        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        out_str: str = ''
        if self.render_note:
            out_str += f'{note_snippet(self.entity)}\n'
        # messages here?
        if self.render_syspath:
            out_str += f'{syspath_snippet()}\n'
        out_str += f'{formatting.format_imports(self.imports)}\n'
        out_str += f'{metadata_snippet(self.entity)}\n'
        out_str += f'{builder_snippet(self.entity)}\n'
        
        # inputs
        out_str += f"{label('INPUTS')}\n"
        out_str += f"{inputs_yaml_snippet()}\n"
        for winp in self.entity.inputs:
            #if not winp.is_runtime:
            out_str += f'{WorkflowInputText(winp).render()}\n'
        out_str += '\n'
        
        # steps (includes tool steps and subworkflow steps)
        for i, step in enumerate(self.entity.steps):
            out_str += f"{label(f'STEP{i+1}: {step.metadata.label}')}\n"
            out_str += '\n'
            out_str += f'{StepText(step, self.entity).render()}\n'
        
        # outputs
        out_str += f"{label('OUTPUTS')}\n\n"
        for wout in self.entity.outputs:
            out_str += f'{WorkflowOutputText(wout).render()}\n'

        if self.render_translate:
            out_str += f'{translate_snippet()}\n'

        return out_str
    
