
from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Tuple

from janis_core.ingestion.galaxy.internal_model.tool.tool import ITool

from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.internal_model.workflow import InputValue, ConnectionInputValue, StaticInputValue, WorkflowInputInputValue
from janis_core.ingestion.galaxy.gxtool.command.components import Positional
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import factory

from .WorkflowInputText import WorkflowInputText
from ..TextRender import TextRender
from .. import formatting
from .. import ordering



### HELPER METHODS ### 

def note() -> str:
    return"""\"\"\"
FOREWORD ----------
Please see [WEBLINK] for information on RUNTIME VALUES
Please see [WEBLINK] for information on the PRE TASK, TOOL STEP, POST TASK structure
\"\"\"
"""

### HELPER CLASSES ###

@dataclass
class ToolInputLine:
    tag_and_value: str 
    special_label: str
    prefix_label: str
    default_label: str
    datatype_label: str

    @cached_property
    # NOTE: I dont like overriding __len__ or any other dunders
    def length(self) -> int:  
        return len(f'{self.tag_and_value},')

    def render(self, info_indent: int) -> str:
        left = f'{self.tag_and_value},'
        right = f'#'
        if self.prefix_label != '':
            right += f' {self.prefix_label}'
        if self.datatype_label != '':
            right += f' [{self.datatype_label}]'
        if self.special_label != '':
            right += f' ({self.special_label})'
        if self.default_label != '':
            right += f' ({self.default_label})'
        justified = f'\t\t{left:<{info_indent+2}}{right}'
        return justified


class ToolInputLineFactory:
    unknown_count = 0
    
    def __init__(self, tool: ITool):
        self.tool = tool

    def create(self, invalue: InputValue) -> ToolInputLine:
        """
        Steps call a tool, and provide a list of key=value pairs
        to feed data to the tool inputs. 

        This func creates a 'key=value,' line to form part of the step call. 
        Feed the correct [value] to the [key] tool input.

        Note: If we have a value but no known tool input, 
        create an 'unknown' dummy tool input to handle situation. 
        """
        if not invalue.component:
            component = self.create_unknown_tool_component(invalue)
            invalue.component = component
        
        tag_and_value = f'{invalue.input_tag}={invalue.wrapped_value}'
        return ToolInputLine(
            tag_and_value=tag_and_value,
            special_label=self.get_special_label(invalue),
            prefix_label=self.get_prefix_label(invalue),
            default_label=self.get_default_label(invalue),
            datatype_label=self.get_datatype_label(invalue),
        )

    def create_unknown_tool_component(self, invalue: InputValue) -> InputComponent:
        """
        create a dummy 'unknown' tool input. 
        Used in situations where we know a step is supplied this InputValue, 
        but don't know which underlying tool input it maps to. 
        """
        self.unknown_count += 1
        name = f'unknown{self.unknown_count}'
        default = None
        if isinstance(invalue, StaticInputValue):
            name = invalue.str_value
            default = invalue.raw_value
        
        positional = factory.positional(name)
        positional.forced_default = default
        positional.forced_docstring = 'ATTENTION: unknown input. please address.'
        self.tool.add_input(positional)
        return positional

    def get_special_label(self, invalue: InputValue) -> str:
        if isinstance(invalue, WorkflowInputInputValue):
            return 'WORKFLOW INPUT'
            # if invalue.is_runtime:
            #     return 'RUNTIME VALUE'
            # else:
            #     return 'WORKFLOW INPUT'
        if isinstance(invalue, ConnectionInputValue):
            return 'CONNECTION'
        return ''
    
    def get_prefix_label(self, invalue: InputValue) -> str:
        if invalue.component:
            if hasattr(invalue.component, 'prefix'):
                return invalue.component.prefix  # type: ignore
            elif isinstance(invalue, Positional):
                return f'position {invalue.cmd_pos}'
        return ''
    
    def get_default_label(self, invalue: InputValue) -> str:
        if isinstance(invalue, StaticInputValue):
            if invalue.is_default:
                return 'DEFAULT'
        return ''
    
    def get_datatype_label(self, invalue: InputValue) -> str:
        if invalue.component:
            return formatting.format_typestr(invalue.component, fmt='value')
        return ''


### MAIN CLASS ###

class StepText(TextRender):
    def __init__(
        self, 
        entity: WorkflowStep, 
        workflow: Workflow,
        render_note: bool=False,
        render_imports: bool=False,
        render_runtime_inputs: bool=False
    ):
        super().__init__()
        self.entity = entity
        self.workflow = workflow
        self.render_note = render_note
        self.render_imports = render_imports
        self.render_runtime_inputs = render_runtime_inputs
    
    @cached_property
    def lines(self) -> list[ToolInputLine]:
        invalues = self.entity.inputs.all
        # comment following line to show all step input values
        invalues = [x for x in invalues if not isinstance(x, StaticInputValue)]
        invalues = ordering.order_step_inputs(invalues)
        factory = ToolInputLineFactory(self.entity.tool)
        return [factory.create(val) for val in invalues]
    
    @cached_property
    def line_len(self) -> int:
        return max([line.length for line in self.lines])

    @property
    def imports(self) -> list[Tuple[str, str]]:
        imports: list[Tuple[str, str]] = []
        if len(self.get_scatter_input_names()) > 1:
            imports.append(('janis_core', 'ScatterDescription'))
            imports.append(('janis_core', 'ScatterMethods'))
        return imports

    def render(self) -> str:
        out_str: str = ''
        if self.render_note:
            out_str += f'{note()}\n'
        if self.render_imports:
            out_str += f'{formatting.format_imports(self.imports)}\n'
        if self.render_runtime_inputs:
            out_str += f'{self.format_runtime_inputs()}\n'
        
        out_str += f'{self.format_step()}\n'
        return out_str

    def format_runtime_inputs(self) -> str:
        inputs = self.entity.inputs.all
        inputs = [x for x in inputs if isinstance(x, WorkflowInputInputValue)]
        inputs = [x for x in inputs if x.is_runtime]
        winps = [self.workflow.get_input(x.input_uuid) for x in inputs]
        out_str = ''
        for winp in winps:
            out_str += f'{WorkflowInputText(winp).render()}\n'
        return out_str

    def format_step(self) -> str:
        step_tag = self.entity.tag
        tool_tag = self.entity.tool.tag
        scatter_stmt = self.format_scatter()
        
        out_str: str = ''
        out_str += 'w.step(\n'
        out_str += f'\t"{step_tag}",\n'
        out_str += f'\t{tool_tag}(\n'
        for line in self.lines:
            out_str += f'{line.render(self.line_len)}\n'
        out_str += '\t)'
        if scatter_stmt:
            out_str += f',\n{scatter_stmt}'
        out_str += '\n)'

        return out_str

    def format_scatter(self) -> Optional[str]:
        names = self.get_scatter_input_names()
        return self.format_scatter_statement(names)

    def get_scatter_input_names(self) -> list[str]:
        # get scatter inputs names
        unknown_count: int = 0
        targets: list[str] = []
        invalues = self.entity.inputs.all
        invalues = ordering.order_step_inputs(invalues)
        for val in invalues:
            if val.scatter:
                if not val.component:
                    unknown_count += 1
                    targets.append(f'{val.input_tag}{unknown_count}')
                else:
                    targets.append(val.input_tag)
        return targets

    def format_scatter_statement(self, names: list[str]) -> Optional[str]:
        # format text
        if len(names) == 0:
            return None
        elif len(names) == 1:
            return f'\tscatter="{names[0]}"'
        else:
            return f'\tscatter=ScatterDescription(fields={repr(names)}, method=ScatterMethods.dot)'


        

