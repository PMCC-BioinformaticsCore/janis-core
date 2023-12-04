

# from typing import Any, Optional

# # tool entities 
# from janis_core.code.codetool import CodeTool
# from janis_core.tool.commandtool import (
#     CommandTool, 
#     CommandToolBuilder, 
#     ToolArgument, 
#     ToolInput, 
#     ToolOutput,
# )
# # # workflow entities
# # from janis_core.graph.steptaginput import Edge, StepTagInput
# # from janis_core.workflow.workflow import (
# #     Workflow,
# #     WorkflowBuilder,
# #     InputNode,
# #     StepNode,
# #     OutputNode,
# # )

# ### PUBLIC ###

# def get_owner_uuid(entity: Any) -> Optional[str]:
#     if isinstance(entity, CommandTool | CommandToolBuilder | CodeTool):
#         return _get_uuid_tool(entity)
#     elif isinstance(entity, ToolInput):
#         return _get_uuid_tool_inp(entity)
#     elif isinstance(entity, ToolArgument):
#         return _get_uuid_tool_arg(entity)
#     elif isinstance(entity, ToolOutput):
#         return _get_uuid_tool_out(entity)
#     else:
#         raise NotImplementedError(f"get_owner_uuid not implemented for this entity")

# ### PRIVATE ###

# # TOOL ENTITIES

# def _get_uuid_tool(entity: CommandTool | CommandToolBuilder | CodeTool) -> Optional[str]:
#     return entity.uuid

# def _get_uuid_tool_arg(entity: ToolArgument) -> Optional[str]:
#     raise NotImplementedError

# def _get_uuid_tool_inp(entity: ToolInput) -> Optional[str]:
#     raise NotImplementedError

# def _get_uuid_tool_out(entity: ToolOutput) -> Optional[str]:
#     raise NotImplementedError


# # WORKFLOW ENTITIES

# # def _get_uuid_workflow(entity: Workflow | WorkflowBuilder) -> Optional[str]:
# #     raise NotImplementedError

# # def _get_uuid_workflow_inp(entity: InputNode) -> Optional[str]:
# #     raise NotImplementedError

# # def _get_uuid_workflow_step(entity: StepNode) -> Optional[str]:
# #     raise NotImplementedError

# # def _get_uuid_workflow_out(entity: OutputNode) -> Optional[str]:
# #     raise NotImplementedError

# # def _get_uuid_workflow_edge(entity: Edge) -> Optional[str]:
# #     raise NotImplementedError

# # def _get_uuid_workflow_sti(entity: StepTagInput) -> Optional[str]:
# #     raise NotImplementedError