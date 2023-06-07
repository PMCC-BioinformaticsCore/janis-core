
from janis_core import Workflow, WorkflowBuilder
from janis_core import CommandTool, CommandToolBuilder
from janis_core import Tool


def to_builders(entity: Tool) -> Tool:
    if isinstance(entity, Workflow):
        for step in entity.step_nodes.values():
            step.tool = to_builders(step.tool)
        if not isinstance(entity, WorkflowBuilder):
            entity = to_workflow_builder(entity)
    elif isinstance(entity, CommandTool) and not isinstance(entity, CommandToolBuilder):
        entity = to_commandtool_builder(entity)
    return entity

def to_workflow_builder(workflow: Workflow) -> WorkflowBuilder:
    # init WorkflowBuilder
    builder = WorkflowBuilder(
        identifier=workflow.id(),
        friendly_name=workflow.friendly_name(),
        version=workflow.version(),
        metadata=workflow.metadata,
        tool_provider=workflow.tool_provider(),
        tool_module=workflow.tool_module(),
        doc=workflow.doc()
    )

    # Add Workflow attributes
    builder.nodes = workflow.nodes
    builder.input_nodes = workflow.input_nodes
    builder.step_nodes = workflow.step_nodes
    builder.output_nodes = workflow.output_nodes
    builder.has_scatter = workflow.has_scatter
    builder.has_subworkflow = workflow.has_subworkflow
    builder.has_multiple_inputs = workflow.has_multiple_inputs

    # Add Tool attributes
    builder.uuid = workflow.uuid
    builder.connections = workflow.connections
    return builder

def to_commandtool_builder(tool: CommandTool) -> CommandToolBuilder:
    return tool.to_command_tool_builder()

