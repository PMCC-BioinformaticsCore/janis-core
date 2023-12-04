
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from janis_core import WorkflowBuilder, CommandToolBuilder, CodeTool
    from janis_core.workflow.workflow import StepNode
    from typing import Any
 
from janis_core.tool.tool import ToolType
from .enums import FormatCategory


def gather_uuids(entity: WorkflowBuilder | CommandToolBuilder | CodeTool | StepNode) -> dict[str, FormatCategory]:
    """
    For this workflow / tool to be rendered, collect the relevant entity uuids
    and assign to categories. 

    This function enables us to use the messages.log_message(entity.uuid) function 
    from anywhere** in the codebase. 
    
    When we are ready to render the workflow / tool, we:
    - Gather all relevant uuids for this file
    - Assign to categories (Tool, ToolInput, Workflow, WorkflowStep etc)
    - Pull messages from the logfile for these uuids 
      [ messages.load_loglines() ]
    - Format and inject the messages in the right spot in the file
      [ inject_messages_tool() or inject_messages_workflow() ]
    
    **except anything with CodeTool. If you want to use messages.log_message() for a 
    CodeTool, 
    for the following entities:
    - CommandTool 
    - CommandToolBuilder 
    - ToolArgument 
    - ToolInput 
    - ToolOutput
    - Workflow
    - WorkflowBuilder
    - InputNode
    - StepNode
    - OutputNode
    - Edge
    - StepTagInput 

    **NOTE: If you want to use messages.log_message() for a CodeTool / PythonTool, you must specify 
    entity_uuid as the CodeTool.uuid rather than a subentity. eg:

        tool = MyPythonTool(PythonTool)
        my_input = tool.inputs()[0]

        messages.log_message(tool.uuid, msg='{my_input.id()} is weird')
                             ^^^^^^^^^ THIS
        
        messages.log_message(my_input.uuid, msg='{my_input.id()} is weird')
                             ^^^^^^^^^^^^^ NOT THIS
    """
    tracer = EntityUuidTracer()
    tracer.trace(entity)
    return tracer.uuids_categories


class EntityUuidTracer:
    """
    Gathers entity uuids for a tool or workflow
    """
    def __init__(self) -> None:
        self.uuids_categories: dict[str, FormatCategory] = {}
    
    def trace(self, entity: Any) -> None:
        if entity.__class__.__name__ == 'StepNode':
            self.trace_step(entity)
        elif entity.type() == ToolType.Workflow:
            self.trace_workflow(entity)
        elif entity.type() == ToolType.CommandTool:
            self.trace_tool(entity)
        elif entity.type() == ToolType.CodeTool:
            self.trace_codetool(entity)
        else:
            raise NotImplementedError

    def trace_tool(self, entity: CommandToolBuilder) -> None:
        self.uuids_categories[entity.uuid] = FormatCategory.MAIN
        
        inputs = entity.inputs() or []
        args = entity.arguments() or []
        outputs = entity.outputs() or []
        
        for inp in inputs:
            self.uuids_categories[inp.uuid] = FormatCategory.INPUT
        for arg in args:
            self.uuids_categories[arg.uuid] = FormatCategory.ARGUMENT
        for out in outputs:
            self.uuids_categories[out.uuid] = FormatCategory.OUTPUT
    
    def trace_codetool(self, entity: CodeTool) -> None:
        self.uuids_categories[entity.uuid] = FormatCategory.MAIN
    
    def trace_workflow(self, entity: WorkflowBuilder) -> None:
        self.uuids_categories[entity.uuid] = FormatCategory.MAIN

        inputs = entity.input_nodes.values() or []
        steps = entity.step_nodes.values() or []
        outputs = entity.output_nodes.values() or []

        # inputs
        for inp in inputs:
            self.uuids_categories[inp.uuid] = FormatCategory.INPUT
        
        # steps
        for step in steps:
            self.trace_step(step)

        # outputs
        for out in outputs:
            self.uuids_categories[out.uuid] = FormatCategory.OUTPUT

    def trace_step(self, entity: StepNode) -> None:
        self.uuids_categories[entity.uuid] = FormatCategory.STEP
        for sti in entity.sources.values():
            self.uuids_categories[sti.uuid] = FormatCategory.STEP
            for edge in sti.source_map:
                self.uuids_categories[edge.uuid] = FormatCategory.STEP



